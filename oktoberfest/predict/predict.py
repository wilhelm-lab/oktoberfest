import logging
import re
from math import ceil
from multiprocessing import current_process
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from spectrum_fundamentals.metrics.similarity import SimilarityMetrics
from tqdm.auto import tqdm
from tritonclient.grpc import InferenceServerClient, InferInput, InferRequestedOutput

from ..data.spectra import FragmentType, Spectra

logger = logging.getLogger(__name__)


def grpc_predict(
    library: Spectra,
    url: str,
    intensity_model: str,
    irt_model: str,
    ssl: bool = True,
    alignment: bool = False,
    job_type: str = "",
):
    """
    Use grpc to predict library and add predictions to library.

    :param library: Spectra object with the library
    :param url: Url including the port of the prediction server
    :param intensity_model: the name of the intensity model on the server
    :param irt_model: the name of the irt model on the server
    :param ssl: whether or not the server requires an ssl encrypted transportation, default = True
    :param alignment: True if alignment present
    :param job_type: TODO
    :return: grpc predictions if we are trying to generate spectral library
    """
    triton_client = InferenceServerClient(url=url, ssl=ssl)
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
    if "tmt" in intensity_model.lower() or "ptm" in intensity_model.lower():
        intensity_input_data["fragmentation_types"] = (
            library.spectra_data["FRAGMENTATION"].to_numpy().reshape(-1, 1).astype(np.object_),
            "BYTES",
        )

    intensity_predictions = infer_predictions(
        triton_client,
        model=intensity_model,
        input_data=intensity_input_data,
        outputs=intensity_outputs,
        batch_size=batch_size,
    )
    intensity_predictions["intensities"][np.where(intensity_predictions["intensities"] < 1e-7)] = 0.0

    irt_input_data = {"peptide_sequences": intensity_input_data["peptide_sequences"]}
    irt_outputs = ["irt"]
    irt_predictions = infer_predictions(
        triton_client,
        model=irt_model,
        input_data=irt_input_data,
        outputs=irt_outputs,
        batch_size=batch_size,
    )

    if job_type == "SpectralLibraryGeneration":
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

    n_batches = ceil(num_spec / batch_size)
    process_identity = current_process()._identity
    if len(process_identity) > 0:
        position = process_identity[0]
    else:
        position = 0

    with tqdm(
        total=n_batches,
        position=position,
        desc=f"Inferring predictions for {num_spec} spectra with batch site {batch_size}",
        leave=True,
    ) as progress:
        for i in range(0, n_batches):
            progress.update(1)
            # logger.info(f"Predicting batch {i+1}/{n_batches}.")
            infer_inputs = []
            for input_key, (data, dtype) in input_data.items():
                batch_data = data[i * batch_size : (i + 1) * batch_size]
                infer_input = InferInput(input_key, batch_data.shape, dtype)
                infer_input.set_data_from_numpy(batch_data)
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


def _prepare_alignment_df(library: Spectra, ce_range: Tuple[int, int], group_by_charge: False) -> Spectra:
    """
    Prepare an alignment DataFrame from the given Spectra library.

    This function creates an alignment DataFrame by removing decoy and HCD fragmented spectra
    from the input library, selecting the top 1000 highest-scoring spectra, and repeating the
    DataFrame for each collision energy (CE) in the given range.

    :param library: the library to be propagated
    :param ce_range: the min and max CE to be propagated for alignment in the dataframe
    :param group_by_charge: if true, select the top 1000 spectra independently for each precursor charge
    :return: a library that is modified according to the description above
    """
    alignment_library = Spectra()
    alignment_library.spectra_data = library.spectra_data.copy()

    # Remove decoy and HCD fragmented spectra
    alignment_library.spectra_data = alignment_library.spectra_data[
        (alignment_library.spectra_data["FRAGMENTATION"] == "HCD") & (~alignment_library.spectra_data["REVERSE"])
    ]
    # Select the 1000 highest scoring or all if there are less than 1000
    temp_df = alignment_library.spectra_data.sort_values(by="SCORE", ascending=False)
    if group_by_charge:
        temp_df = temp_df.groupby("PRECURSOR_CHARGE")

    alignment_library.spectra_data = temp_df.head(1000)

    # Repeat dataframe for each CE
    ce_range = range(*ce_range)
    nrow = len(alignment_library.spectra_data)
    alignment_library.spectra_data = pd.concat([alignment_library.spectra_data for _ in ce_range], axis=0)
    alignment_library.spectra_data["ORIG_COLLISION_ENERGY"] = alignment_library.spectra_data["COLLISION_ENERGY"]
    alignment_library.spectra_data["COLLISION_ENERGY"] = np.repeat(ce_range, nrow)
    alignment_library.spectra_data.reset_index(inplace=True)
    return alignment_library


def ce_calibration(library: Spectra, ce_range: Tuple[int, int], group_by_charge: bool, **server_kwargs) -> pd.Series:
    """
    Calculate best collision energy for peptide property predictions.

    The function propagates the provided library object to test NCEs in the given ce range, performs
    intensity prediction for the 1000 highest scoring target PSMs at each NCE and computes the spectral angle
    between predicted and observed intensities before returning the alignment library.

    :param library: spectral library to perform CE calibration on
    :param ce_range: the min and max CE to be tested during calibration
    :param group_by_charge: if true, select the top 1000 spectra independently for each precursor charge
    :param server_kwargs: Additional parameters that are forwarded to grpc_predict
    :return: pandas series containing the spectral angle for all tested collision energies
    """
    alignment_library = _prepare_alignment_df(library, ce_range=ce_range, group_by_charge=group_by_charge)
    grpc_predict(alignment_library, alignment=True, **server_kwargs)
    _alignment(alignment_library)
    return alignment_library


def _alignment(alignment_library: Spectra):
    """
    Perform the alignment of predicted versus raw intensities.

    The function calculates the spectral angle between predicted and observed fragment intensities and
    adds it as a column to the alignment library.

    :param alignment_library: the library to perform the alignment on
    """
    pred_intensity = alignment_library.get_matrix(FragmentType.PRED)
    raw_intensity = alignment_library.get_matrix(FragmentType.RAW)
    # return pred_intensity.toarray(), raw_intensity.toarray()
    sm = SimilarityMetrics(pred_intensity, raw_intensity)
    alignment_library.spectra_data["SPECTRAL_ANGLE"] = sm.spectral_angle(raw_intensity, pred_intensity, 0)
    alignment_library.spectra_data = alignment_library.spectra_data[
        alignment_library.spectra_data["SPECTRAL_ANGLE"] != 0
    ]
