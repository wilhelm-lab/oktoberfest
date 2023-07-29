import logging
import re
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from spectrum_fundamentals.metrics.similarity import SimilarityMetrics
from tritonclient.grpc import InferenceServerClient, InferInput, InferRequestedOutput

from ..data.spectra import FragmentType, Spectra
from ..plotting import plot_mean_sa_ce, plot_violin_sa_ce
from ..preprocessing import merge_mzml_and_msms, prepare_alignment_df
from ..utils.config import Config

logger = logging.getLogger(__name__)


def grpc_predict(config: Config, library: Spectra, alignment: bool = False):
    """
    Use grpc to predict library and add predictions to library.

    :param library: Spectra object with the library
    :param alignment: True if alignment present
    :return: grpc predictions if we are trying to generate spectral library
    """
    triton_client = InferenceServerClient(url=config.prediction_server, ssl=config.ssl)
    batch_size = 1000

    intensity_outputs = ["intensities", "mz", "annotation"]
    intensity_input_data = {
        "peptide_sequences": (
            library.spectra_data["MODIFIED_SEQUENCE"].to_numpy().reshape(-1, 1).astype(np.object_),
            "BYTES",
        ),
        "collision_energies": (
            library.spectra_data["COLLISION_ENERGY"].to_numpy().reshape(-1, 1).astype(np.float32),
            "FP32",
        ),
        "precursor_charges": (
            library.spectra_data["PRECURSOR_CHARGE"].to_numpy().reshape(-1, 1).astype(np.int32),
            "INT32",
        ),
    }
    intensity_model = config.models["intensity"]
    if "tmt" in intensity_model.lower() or "ptm" in intensity_model.lower():
        library.spectra_data["FRAGMENTATION_GRPC"] = library.spectra_data["FRAGMENTATION"].apply(
            lambda x: 2 if x == "HCD" else 1
        )
        intensity_input_data["fragmentation_types"] = (
            library.spectra_data["FRAGMENTATION_GRPC"].to_numpy().reshape(-1, 1).astype(np.float32),
            "FP32",
        )

    intensity_predictions = infer_predictions(
        triton_client,
        model=intensity_model,
        input_data=intensity_input_data,
        outputs=intensity_outputs,
        batch_size=batch_size,
    )
    intensity_predictions["intensities"][np.where(intensity_predictions["intensities"] < 1e-7)] = 0.0

    irt_model = config.models["irt"]
    irt_input_data = {"peptide_sequences": intensity_input_data["peptide_sequences"]}
    irt_outputs = ["irt"]
    irt_predictions = infer_predictions(
        triton_client,
        model=irt_model,
        input_data=irt_input_data,
        outputs=irt_outputs,
        batch_size=batch_size,
    )

    if config.job_type == "SpectralLibraryGeneration":
        intensity_prediction_dict = {
            "intensity": intensity_predictions["intensities"],
            "fragmentmz": intensity_predictions["mz"],
            "annotation": parse_fragment_labels(
                intensity_predictions["annotation"],
                library.spectra_data["PRECURSOR_CHARGE"].to_numpy()[:, None],
                library.spectra_data["PEPTIDE_LENGTH"].to_numpy()[:, None],
            ),
        }
        output_dict = {intensity_model: intensity_prediction_dict, irt_model: irt_predictions["irt"]}
        return output_dict

    intensities_pred = pd.DataFrame()
    intensities_pred["intensity"] = intensity_predictions["intensities"].tolist()
    library.add_matrix(intensities_pred["intensity"], FragmentType.PRED)

    if alignment:
        return

    library.add_column(irt_predictions["irt"], name="PREDICTED_IRT")


def infer_predictions(
    triton_client: InferenceServerClient,
    model: str,
    input_data: Dict[str, Tuple[np.ndarray, str]],
    outputs: List[str],
    batch_size: int,
) -> Dict[str, np.ndarray]:
    """
    Infer predictions from a triton client.

    :param triton_client: An inference client using grpc
    :param model: a model that is recognized by the server specified in the triton client
    :param input_data: a dictionary that contains the input names (key) for the specific model
        and a tuple of the input_data as a numpy array of shape [:, 1] and the dtype recognized
        by the triton client (value).
    :param outputs: a list of output names for the specific model
    :param batch_size: the number of elements from the input_data that should be provided to the
        triton client at once
    :return: a dictionary containing the predictions (values) for the given outputs (keys)
    """
    num_spec = len(input_data[list(input_data)[0]][0])
    predictions: Dict[str, List[np.ndarray]] = {output: [] for output in outputs}

    for i in range(0, num_spec, batch_size):
        if num_spec < i + batch_size:
            current_batchsize = num_spec - i
        else:
            current_batchsize = batch_size

        infer_inputs = []
        for input_key, (data, dtype) in input_data.items():
            infer_input = InferInput(input_key, [current_batchsize, 1], dtype)
            infer_input.set_data_from_numpy(data[i : i + current_batchsize])
            infer_inputs.append(infer_input)

        infer_outputs = [InferRequestedOutput(output) for output in outputs]

        prediction = triton_client.infer(model, inputs=infer_inputs, outputs=infer_outputs)

        for output in outputs:
            predictions[output].append(prediction.as_numpy(output))

    return {key: np.vstack(value) for key, value in predictions.items()}


def parse_fragment_labels(spectra_labels: np.ndarray, precursor_charges: np.ndarray, seq_lengths: np.ndarray):
    """Uses regex to parse labels."""
    pattern = rb"([y|b])([0-9]{1,2})\+([1-3])"
    fragment_types = []
    fragment_numbers = []
    fragment_charges = []
    for spectrum_labels in spectra_labels:
        types = []
        numbers = []
        charges = []
        for label in spectrum_labels:
            match = re.match(pattern, label)
            if match:
                groups = match.groups()
                types.append(groups[0].decode())
                numbers.append(int(groups[1]))
                charges.append(int(groups[2]))
            else:
                raise ValueError(f"String {label} does not match the expected fragment label pattern")
        fragment_types.append(types)
        fragment_numbers.append(numbers)
        fragment_charges.append(charges)

    fragment_type_array = np.array(fragment_types)
    fragment_number_array = np.array(fragment_numbers)
    fragment_charge_array = np.array(fragment_charges)
    mask = np.where((fragment_charge_array > precursor_charges) | (fragment_number_array >= seq_lengths))
    fragment_type_array[mask] = "N"
    fragment_number_array[mask] = 0
    fragment_charge_array[mask] = 0

    return {"type": fragment_type_array, "number": fragment_number_array, "charge": fragment_charge_array}


def perform_alignment(config: Config, library: Spectra, df_search: pd.DataFrame) -> int:
    """
    Perform alignment and get the best CE.

    :param df_search: search result as pd.DataFrame
    """
    hdf5_path = get_hdf5_path(config)
    logger.info(f"Path to hdf5 file with annotations for {config.output}: {hdf5_path}")
    if hdf5_path.is_file():
        library.read_from_hdf5(hdf5_path)
    else:
        library = merge_mzml_and_msms(config, library, df_search)
        library.write_as_hdf5(hdf5_path)  # write_metadata_annotation
    if (library.spectra_data["FRAGMENTATION"] == "HCD").any():
        alignment_library = prepare_alignment_df(library)
        grpc_predict(config, alignment_library, alignment=True)  # predict alignment
        ce_alignment = _alignment(config, alignment_library)
        best_ce = ce_alignment.idxmax()
    else:
        best_ce = 35
    return best_ce


def get_hdf5_path(config: Config) -> Path:
    """Get path to hdf5 file."""
    return config.output / "data" / config.spectra.with_suffix(".mzML.hdf5").name


def _alignment(config: Config, alignment_library: Spectra):
    """
    Edit library to try different ranges of ce 15-50. then predict with the new library.

    Check https://gitlab.lrz.de/proteomics/prosit_tools/oktoberfest/-/blob/develop/oktoberfest/ce_calibration/grpc_alignment.py
    """
    pred_intensity = alignment_library.get_matrix(FragmentType.PRED)
    raw_intensity = alignment_library.get_matrix(FragmentType.RAW)
    # return pred_intensity.toarray(), raw_intensity.toarray()
    sm = SimilarityMetrics(pred_intensity, raw_intensity)
    alignment_library.spectra_data["SPECTRAL_ANGLE"] = sm.spectral_angle(raw_intensity, pred_intensity, 0)
    alignment_library.spectra_data.to_csv(config.output / "results" / "SA.tsv", sep="\t")
    ce_alignment = alignment_library.spectra_data.groupby(by=["COLLISION_ENERGY"])["SPECTRAL_ANGLE"].mean()

    plot_mean_sa_ce(
        sa_ce_df=ce_alignment,
        filename=config.output / "results" / f"{config.spectra.stem}_mean_spectral_angle_ce.svg",
        best_ce=ce_alignment.idxmax(),
    )
    plot_violin_sa_ce(
        df=alignment_library.spectra_data[["COLLISION_ENERGY", "SPECTRAL_ANGLE"]],
        filename=config.output / "results" / f"{config.spectra.stem}_violin_spectral_angle_ce.svg",
        best_ce=ce_alignment.idxmax(),
    )
    return ce_alignment
