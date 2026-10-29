import datetime
import json
import logging
import pickle
import shutil
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from functools import partial
from math import ceil
from multiprocessing import Manager, Process, pool
from pathlib import Path
from typing import Optional, Union

import anndata
import numpy as np
import pandas as pd
import scipy.sparse
import spectrum_fundamentals.constants as c
from sklearn.linear_model import LinearRegression, RANSACRegressor
from spectrum_io.spectral_library import MSP, DLib, SpectralLibrary, Spectronaut
from tqdm.auto import tqdm

from oktoberfest import __copyright__, __version__
from oktoberfest import plotting as pl
from oktoberfest import predict as pr
from oktoberfest import preprocessing as pp
from oktoberfest import rescore as re

from .data.spectra import Spectra
from .utils import Config, JobPool, ProcessStep, apply_quant, group_iterator
from .utils.config import ADJACENT_DECOY_RESIDUES

logger = logging.getLogger(__name__)

# column of a library from previous runs (reuseLibraries): collision energy its intensities were predicted at
REUSED_CE_COLUMN = "REUSED_CE"


def _make_predictions_error_callback(failure_progress_tracker, failure_lock, error):
    logger.error(
        f"Prediction failed due to: {error} Batch will be missing from output file. "
        "DO NOT STOP THIS RUN: The index of the batch is stored and your output file will be appended "
        "by the missing batches if you rerun without changing your config file after this run is completed."
    )
    with failure_lock:
        failure_progress_tracker.value += 1


def _make_predictions(config, queue_out, progress, lock, batch_df):
    predictors = {model_key: pr.Predictor.from_config(config, model_type=model_key) for model_key in config.models}
    predictions = {
        output_name: output
        for predictor in predictors.values()
        for output_name, output in predictor._predictor.predict(batch_df).items()
    }
    queue_out.put((predictions, batch_df))
    with lock:
        progress.value += 1


def _preprocess(spectra_files: list[Path], config: Config) -> list[Path]:
    preprocess_search_step = ProcessStep(config.output, "preprocessing_search")
    if not preprocess_search_step.is_done():
        # load search results
        logger.info(f"Converting search results from {config.search_results} to internal search result.")
        msms_output = config.output / "msms"
        msms_output.mkdir(exist_ok=True)
        internal_search_file = msms_output / "msms.prosit"
        tmt_label = config.tag
        # PTM site permutations are only generated if localization is requested, otherwise PTMs are kept as
        # localized by the search engine
        # including decoy residues used to estimate the false localization rate
        ptm_unimods = config.ptm_permutation_sites
        ptm_sites = None
        ptm_model = config.ptm_model
        search_inputs = config.search_inputs
        converted = [
            pp.convert_search(
                input_path=search_results_path,
                search_engine=search_engine,
                tmt_label=tmt_label,
                custom_mods=config.custom_to_unimod(),
                # several search engines: one internal file per engine, merged into msms.prosit
                output_file=(
                    internal_search_file
                    if len(search_inputs) == 1
                    else msms_output / f"msms_{i}_{search_engine}.prosit"
                ),
                ptm_unimod_id=ptm_unimods,
                ptm_sites=ptm_sites,
                ptm_model=ptm_model,
                num_threads=config.num_threads,
            )
            for i, (search_results_path, search_engine) in enumerate(search_inputs)
        ]
        if len(search_inputs) == 1:
            search_results = converted[0]
        else:
            search_results = pp.merge_search_results(
                converted,
                [search_engine for _, search_engine in search_inputs],
                engines_file=msms_output / "search_engines.csv",
            )
            search_results.to_csv(internal_search_file, index=False)
        if config.spectra_type.lower() in ["d", "hdf"]:
            timstof_metadata = pp.convert_timstof_metadata(
                input_path=config.search_results,
                search_engine=config.search_results_type,
                output_file=msms_output / "tims_meta.csv",
            )

        # TODO add support for internal timstof metadata
        logger.info(f"Read {len(search_results)} PSMs from {internal_search_file}")
        if config.ptm_decoy_residues and config.ptm_decoy_selection == "adjacent":
            search_results = re.add_adjacent_decoy_candidates(
                search_results, config.ptm_sites, config.ptm_decoy_priority_residues, ADJACENT_DECOY_RESIDUES
            )
        elif config.ptm_decoy_residues:
            search_results = re.select_decoy_candidates(search_results, config.ptm_sites, config.ptm_decoy_residues)
        if config.ptm_decoy_residues:
            if config.ptm_localization_options.get("decoy_proxy") == "search":
                fractions = re.decoy_proxy_fractions(search_results, config.ptm_sites)
                config.decoy_proxy_fractions_file.write_text(json.dumps(fractions))

        # Filter peptides only if using Koina intensity model
        if not config.models.get("dlomix_intensity"):
            model_type = config.models.get("intensity", "")
            try:
                search_results = pp.filter_peptides_for_model(peptides=search_results, model=model_type)
            except ValueError:
                logger.exception(
                    ValueError(
                        f"Unknown model {model_type}. Please ensure it is one of ['prosit', 'ms2pip', 'alphapept']."
                        "If you're using local prediction, please ensure the model type is contained in the model file name."
                    )
                )

        # split search results
        searchfiles_found = pp.split_search(
            search_results=search_results,
            output_dir=config.output / "msms",
            filenames=[spectra_file.stem for spectra_file in spectra_files],
        )
        # split timstof metadata
        if config.spectra_type.lower() in ["d", "hdf"]:
            _ = pp.split_timstof_metadata(
                timstof_metadata=timstof_metadata,
                output_dir=config.output / "msms",
                filenames=searchfiles_found,
            )
        preprocess_search_step.mark_done()
    else:
        searchfiles_found = [msms_file.stem for msms_file in (config.output / "msms").glob("*rescore")]
    spectra_files_to_return = []
    for spectra_file in spectra_files:
        if spectra_file.stem in searchfiles_found:
            spectra_files_to_return.append(spectra_file)

    return spectra_files_to_return


def _get_annotation_fragmentation(spectra: pd.DataFrame, config: Config) -> str:
    """
    Get the fragmentation method used to annotate the fragment ions of a spectra file.

    The fragmentation_method in the config takes precedence. Otherwise, the fragmentation read from the spectra file
    is used if all spectra share one supported method. Unknown, unsupported or mixed fragmentation methods fall back
    to HCD with a warning.

    :param spectra: the spectra read from the spectra file
    :param config: the configuration object
    :return: the fragmentation method used for annotation
    """
    if config.data.get("fragmentation_method") is not None:
        return config.fragmentation_method
    # TODO: refine annotation for files with unknown, unsupported or mixed fragmentation methods, e.g. by annotating
    #  each method with its own ion types instead of falling back to HCD
    methods = spectra["FRAGMENTATION"].astype(str).str.upper().unique()
    if len(methods) == 1 and methods[0] in c.FRAGMENTATION_TO_IONS_BY_PAIRS:
        return methods[0]
    logger.warning(
        f"Found fragmentation methods {sorted(methods)} in the spectra, which are unknown, unsupported or mixed. "
        "Fragment ions are annotated as HCD. Set fragmentation_method in the config to annotate them differently."
    )
    return "HCD"


def _annotate_and_get_library(spectra_file: Path, config: Config, tims_meta_file: Optional[Path] = None) -> Spectra:
    data_dir = config.output / "data"
    data_dir.mkdir(exist_ok=True)
    hdf5_path = data_dir / spectra_file.with_suffix(".mzml.hdf5").name
    if hdf5_path.is_file():
        aspec = Spectra.from_hdf5(hdf5_path)
        instrument_type = config.instrument_type
        if instrument_type is not None and aspec.obs["INSTRUMENT_TYPES"].values[0] != instrument_type:
            aspec.obs["INSTRUMENT_TYPES"] = instrument_type
            aspec.write_as_hdf5(hdf5_path)
    else:
        spectra_dir = config.output / "spectra"
        spectra_dir.mkdir(exist_ok=True)
        format_ = spectra_file.suffix.lower()
        if format_ == ".raw":
            file_to_load = spectra_dir / spectra_file.with_suffix(".mzML").name
            pp.convert_raw_to_mzml(spectra_file, file_to_load, thermo_exe=config.thermo_exe)
        elif format_ in [".mzml", ".hdf"]:
            file_to_load = spectra_file
        elif format_ == ".d":
            file_to_load = spectra_dir / spectra_file.with_suffix(".hdf").name
            pp.convert_d_to_hdf(spectra_file, file_to_load)
        spectra = pp.load_spectra(file_to_load, tims_meta_file=tims_meta_file)
        config_instrument_type = config.instrument_type
        if config_instrument_type is not None:
            spectra["INSTRUMENT_TYPES"] = config_instrument_type
        # only override the per-spectrum fragmentation read from the spectra file if explicitly configured
        if config.data.get("fragmentation_method") is not None:
            spectra["FRAGMENTATION"] = config.fragmentation_method
        search = pp.load_search(config.output / "msms" / spectra_file.with_suffix(".rescore").name)
        library = pp.merge_spectra_and_peptides(spectra, search)

        # Check for crosslinked or local model (skip for DLOmix)
        intensity_model = config.models.get("intensity", "").lower()
        if "xl" in intensity_model:
            if "cms2" in intensity_model:
                cms2 = True
            else:
                cms2 = False
            aspec = pp.annotate_spectral_library_xl(
                psms=library, cms2=cms2, mass_tol=config.mass_tolerance, unit_mass_tol=config.unit_mass_tolerance
            )
            if config.search_results_type == "xisearch":
                aspec.obs["start_pos_p1"] = aspec.obs["start_pos_p1"].astype(str)
                aspec.obs["start_pos_p2"] = aspec.obs["start_pos_p2"].astype(str)
        else:
            annotate_neutral_loss = config.ptm_use_neutral_loss
            aspec = pp.annotate_spectral_library(
                psms=library,
                mass_tol=config.mass_tolerance,
                unit_mass_tol=config.unit_mass_tolerance,
                fragmentation_method=_get_annotation_fragmentation(library, config),
                p_window=config.p_window,
                custom_mods=config.unimod_to_mass(),
                annotate_neutral_loss=annotate_neutral_loss,
                multifrag=config.check_multifrag(),
                featured_ions=config.ion_types,
                phospho_loss_residues=config.ptm_phospho_loss_residues,
                plain_loss_residues=config.ptm_plain_loss_residues,
            )

        aspec.write_as_hdf5(hdf5_path)  # write_metadata_annotation

    return aspec


def _ce_calibration_library(library: Spectra, config: Config) -> Spectra:
    """
    Get the PSMs used for CE calibration.

    If PTMs are localized, the library contains all site permutations of each PSM. They share the search engine score
    of the original PSM, so the PSMs selected for calibration would mostly be wrong localizations. Only PSMs without
    the localized PTMs are used for calibration in this case.

    :param library: the library with all PSMs
    :param config: the configuration object
    :return: the PSMs to use for CE calibration
    """
    if not config.ptm_sites:
        return library
    ptm_pattern = r"UNIMOD:(?:" + "|".join(str(ptm) for ptm in config.ptm_sites) + r")\]"
    without_ptm = ~library.obs["MODIFIED_SEQUENCE"].str.contains(ptm_pattern, regex=True).to_numpy(dtype=bool)
    if not without_ptm.any():
        logger.warning(
            "No PSMs without the localized PTMs found, using all PSMs including permutations for CE calibration."
        )
        return library
    logger.info(f"Using {without_ptm.sum()} PSMs without the localized PTMs for CE calibration.")
    return Spectra(library[without_ptm].copy())


def _map_best_ce(library: Spectra, best_ce_per_orig: pd.Series, is_hcd: np.ndarray, library_column: str):
    # map from a snapshot of the original values in one step, so a CE that was just assigned is never
    # matched again by another group (e.g. 30 -> 33 must not then be remapped by the 33 -> 35 group).
    # Only HCD spectra are calibrated, so CID spectra set to 35 above are never remapped, even if an
    # HCD group also had an original CE of 35. HCD groups without a calibration keep their original CE.
    orig_values = library.obs.loc[is_hcd, library_column]
    calibrated_ce = orig_values.map(best_ce_per_orig)
    n_uncalibrated = int(calibrated_ce.isna().sum())
    if n_uncalibrated > 0:
        logger.warning(
            f"{n_uncalibrated} HCD spectra have no calibrated CE for their {library_column} and keep their original CE"
        )
    library.obs.loc[is_hcd, "COLLISION_ENERGY"] = calibrated_ce.fillna(
        library.obs.loc[is_hcd, "COLLISION_ENERGY"]
    ).to_numpy()


def _stored_best_ce(spectra_file: Path, config: Config) -> Optional[pd.Series]:
    """
    Get the best CE per original CE (DDA) or per precursor charge (DIA) from a previous calibration.

    :param spectra_file: the spectra file the calibration belongs to
    :param config: the configuration object; ce_alignment_options.ce_calibration_dir points to the previous results
    :return: the best CE per group, or None if no previous calibration is configured or found for this file
    """
    ce_calibration_dir = config.ce_calibration_dir
    if ce_calibration_dir is None:
        return None
    ce_file = ce_calibration_dir / f"{spectra_file.stem}_ce.txt"
    if not ce_file.is_file():
        logger.warning(f"No previous CE calibration found for {spectra_file.stem} in {ce_calibration_dir}")
        return None
    with open(ce_file) as f:
        best_ce = json.load(f)
    if not isinstance(best_ce, dict):
        return None
    to_key = int if config.acquisition_method.lower() == "dia" else float
    return pd.Series({to_key(float(group)): ce for group, ce in best_ce.items()})


def _get_best_ce(library: Spectra, spectra_file: Path, config: Config):
    results_dir = config.output / "results/ce_calibration"
    results_dir.mkdir(parents=True, exist_ok=True)
    # the calibrated CE only applies to HCD spectra; all other spectra (e.g. CID) are set to 35
    # same selection as in the CE alignment, which only uses HCD spectra
    is_hcd = (library.obs["FRAGMENTATION"] == "HCD").to_numpy(dtype=bool)
    library.obs.loc[~is_hcd, "COLLISION_ENERGY"] = 35
    stored_best_ce = None
    if is_hcd.any() and not config.use_ransac_model:
        stored_best_ce = _stored_best_ce(spectra_file, config)
    if stored_best_ce is not None:
        logger.info(f"Using the previous CE calibration for {spectra_file.stem}: {stored_best_ce.to_dict()}")
        library_column = "PRECURSOR_CHARGE" if config.acquisition_method.lower() == "dia" else "COLLISION_ENERGY"
        _map_best_ce(library, stored_best_ce, is_hcd, library_column)
        with open(results_dir / f"{spectra_file.stem}_ce.txt", "w") as f:
            json.dump(stored_best_ce.to_dict(), f)
    elif is_hcd.any():
        use_ransac_model = config.use_ransac_model
        if config.acquisition_method.lower() == "dia" or use_ransac_model:
            group_by_charge = True
        else:
            group_by_charge = False
        predictor = pr.Predictor.from_config(config, model_type="intensity")

        # Check for crosslinked model (skip for DLOmix)
        intensity_model = config.models.get("intensity", "").lower()
        is_xl = "xl" in intensity_model

        if is_xl:
            alignment_library = predictor.ce_calibration(
                library,
                config.ce_range,
                use_ransac_model,
                xl=True,
                dataset_name=spectra_file.stem + "_ce_calibration",
            )
        else:
            alignment_library = predictor.ce_calibration(
                _ce_calibration_library(library, config),
                config.ce_range,
                group_by_charge,
                dataset_name=spectra_file.stem + "_ce_calibration",
            )

        if use_ransac_model:
            logger.info("Performing RANSAC regression")
            calib_group = (
                alignment_library.obs.groupby(
                    by=["PRECURSOR_CHARGE", "ORIG_COLLISION_ENERGY", "COLLISION_ENERGY", "MASS"], as_index=False
                )["SPECTRAL_ANGLE"]
                .mean()
                .groupby(["PRECURSOR_CHARGE", "ORIG_COLLISION_ENERGY", "MASS"], as_index=False)
                .apply(lambda x: x.loc[x["SPECTRAL_ANGLE"].idxmax()])
            )
            calib_group["delta_collision_energy"] = (
                calib_group["COLLISION_ENERGY"] - calib_group["ORIG_COLLISION_ENERGY"]
            )
            x = calib_group[["MASS", "PRECURSOR_CHARGE"]]  # input feature
            y = calib_group["delta_collision_energy"]  # target variable
            ransac = RANSACRegressor(LinearRegression(), residual_threshold=1.5, random_state=42)
            ransac.fit(x, y)

            for charge, df in calib_group.groupby("PRECURSOR_CHARGE"):
                r2_score = ransac.score(df[["MASS", "PRECURSOR_CHARGE"]], df["COLLISION_ENERGY"])
                title = f"Scatter Plot with RANSAC Model \nSlope: {ransac.estimator_.coef_[0]:.2f}, "
                title += f"Intercept: {ransac.estimator_.intercept_:.2f}, R2: {r2_score:.2f}"
                pl.plot_ce_ransac_model(
                    sa_ce_df=df,
                    filename=results_dir / f"{spectra_file.stem}_ce_ransac_model_{charge}.svg",
                    title=title,
                )

            delta_ce = ransac.predict(library.obs.loc[is_hcd, ["MASS", "PRECURSOR_CHARGE"]])
            library.obs.loc[is_hcd, "COLLISION_ENERGY"] = np.maximum(
                0, library.obs.loc[is_hcd, "COLLISION_ENERGY"] + delta_ce
            )

        else:
            # DIA: one calibrated CE per precursor charge; DDA: one calibrated CE per original CE.
            # library_column is the column in library.obs holding the same grouping values; for DDA the
            # library's COLLISION_ENERGY still holds the original CE (ORIG_COLLISION_ENERGY only exists
            # in alignment_library)
            if config.acquisition_method.lower() == "dia":
                groupby_column = "PRECURSOR_CHARGE"
                library_column = "PRECURSOR_CHARGE"
            else:
                groupby_column = "ORIG_COLLISION_ENERGY"
                library_column = "COLLISION_ENERGY"
            ce_alignment = alignment_library.obs.groupby(by=[groupby_column, "COLLISION_ENERGY"])[
                "SPECTRAL_ANGLE"
            ].mean()

            best_ce_per_orig = ce_alignment.groupby(level=groupby_column).idxmax().apply(lambda x: x[1])
            logger.info(f"Best CE per original CE: {best_ce_per_orig.to_dict()}")

            for orig_ce in ce_alignment.groupby(level=groupby_column).groups.keys():
                pl.plot_mean_sa_ce(
                    sa_ce_df=ce_alignment.to_frame().reset_index()[
                        ce_alignment.to_frame().reset_index()[groupby_column] == orig_ce
                    ],
                    filename=results_dir / f"{spectra_file.stem}_mean_spectral_angle_ce_{orig_ce}.svg",
                )
                pl.plot_violin_sa_ce(
                    sa_ce_df=alignment_library.obs[alignment_library.obs[groupby_column] == orig_ce][
                        ["COLLISION_ENERGY", "SPECTRAL_ANGLE"]
                    ],
                    filename=results_dir / f"{spectra_file.stem}_violin_spectral_angle_ce_{orig_ce}.svg",
                )

            _map_best_ce(library, best_ce_per_orig, is_hcd, library_column)
            with open(results_dir / f"{spectra_file.stem}_ce.txt", "w") as f:
                json.dump(best_ce_per_orig.to_dict(), f)
                f.close()

    else:
        best_ce = 35
        library.obs["COLLISION_ENERGY"] = best_ce

        with open(results_dir / f"{spectra_file.stem}_ce.txt", "w") as f:
            json.dump(best_ce, f)
            f.close()


def _speclib_from_digestion(config: Config) -> Spectra:
    library_input_type = config.library_input_type
    peptide_dict = None
    library_file = config.library_input
    if library_input_type == "fasta":
        p2p_file = config.output / "peptide_to_proteins.csv"
        digest_step = ProcessStep(config.output, "speclib_digested")
        if not digest_step.is_done():
            peptide_dict = pp.digest(
                fasta=config.library_input,
                digestion=config.digestion,
                missed_cleavages=config.missed_cleavages,
                db=config.db,
                enzyme=config.enzyme,
                special_aas=config.special_aas,
                min_length=config.min_length,
                max_length=config.max_length,
            )
            # Convert dictionary to DataFrame
            p2p_df = pd.DataFrame(list(peptide_dict.items()), columns=["peptide", "proteins"])
            p2p_df["proteins"] = p2p_df["proteins"].apply(lambda x: ";".join(x))
            p2p_df.to_csv(p2p_file, index=False)
            digest_step.mark_done()
        library_input_type = "peptides"
        library_file = p2p_file

    if library_input_type == "peptides":
        internal_library_file = config.output / "peptides_internal.csv"
        created_internal_step = ProcessStep(config.output, "speclib_created_internal")
        if not created_internal_step.is_done():
            proteins = None

            if peptide_dict is None:
                p2p_df = pd.read_csv(library_file)
                if "proteins" in p2p_df.columns:
                    p2p_df["proteins"].fillna("unknown", inplace=True)
                    proteins = p2p_df["proteins"].apply(lambda x: x.split(";")).to_list()
                peptides = p2p_df["peptide"].to_list()
            else:
                peptides = list(peptide_dict.keys())
                proteins = list(peptide_dict.values())
            internal_df = pp.generate_metadata(
                peptides=peptides,
                collision_energy=config.collision_energy,
                precursor_charge=config.precursor_charge,
                fragmentation=config.fragmentation,
                nr_ox=config.nr_ox,
                instrument_type=config.instrument_type,
                proteins=proteins,
            )
            library_file = config.output / "peptides_internal.csv"
            internal_df.to_csv(internal_library_file, sep=",", index=False)
            created_internal_step.mark_done()
        library_file = internal_library_file

    else:
        raise ValueError(f'Library input type {library_input_type} not understood. Can only be "fasta", "peptides".')
    spec_library = pp.gen_lib(library_file)

    pp_and_filter_step = ProcessStep(config.output, "speclib_filtered")

    data_dir = config.output / "data"
    if not pp_and_filter_step.is_done():
        data_dir.mkdir(exist_ok=True)
        spec_library = pp.process_and_filter_spectra_data(
            library=spec_library, model=config.models.get("intensity", ""), tmt_label=config.tag
        )
        spec_library.write_as_hdf5(data_dir / f"{library_file.stem}_filtered.hdf5")
        pp_and_filter_step.mark_done()
    else:
        spec_library = Spectra.from_hdf5(data_dir / f"{library_file.stem}_filtered.hdf5")

    return spec_library


def _get_writer_and_output(results_path: Path, output_format: str) -> tuple[type[SpectralLibrary], Path]:
    libfile_prefix = "predicted_library"
    if output_format == "msp":
        return MSP, results_path / f"{libfile_prefix}.msp"
    elif output_format == "spectronaut":
        return Spectronaut, results_path / f"{libfile_prefix}.csv"
    elif output_format == "dlib":
        return DLib, results_path / f"{libfile_prefix}.dlib"
    else:
        raise ValueError(f"{output_format} is not supported as spectral library type")


def _get_batches_and_mode(out_file: Path, failed_batch_file: Path, obs: pd.DataFrame, batch_size: int, model: str):
    if out_file.is_file():
        if failed_batch_file.is_file():
            with open(failed_batch_file, "rb") as fh:
                batch_iterator = pickle.load(fh)
            mode = "a"
            logger.warning(
                f"Found existing spectral library {out_file}. Attempting to append missing batches from previous run..."
            )
        else:
            logger.error(
                f"A file {out_file} already exists but no information about missing batches "
                "from a previous run could be found. Stopping to prevent corruption / data loss. "
                "If this is intended, delete the file and rerun."
            )
            sys.exit(1)
    else:
        if "alphapept" in model.lower():
            batch_iterator = group_iterator(df=obs, group_by_column="PEPTIDE_LENGTH", max_batch_size=batch_size)
        else:
            batch_iterator = (
                obs.index[i * batch_size : (i + 1) * batch_size].to_numpy() for i in range(ceil(len(obs) / batch_size))
            )
        mode = "w"

    return list(batch_iterator), mode


def _update(pbar: tqdm, postfix_values: dict[str, int]):
    total_val = sum(postfix_values.values())
    if total_val > pbar.n:
        pbar.set_postfix(**postfix_values)
        pbar.n = total_val
        pbar.refresh()


def _check_write_failed_batch_file(failed_batch_file: Path, n_failed: int, results: list[pool.AsyncResult]):
    if n_failed > 0:
        failed_batches = []
        for i, result in enumerate(results):
            try:
                result.get()
            except Exception:
                failed_batches.append(i)
        logger.error(
            f"Prediction for {n_failed} / {i + 1} batches failed. Check the log to find out why. "
            "Then rerun without changing the config file to append only the missing batches to your output file."
        )
        with open(failed_batch_file, "wb") as fh:
            pickle.dump(failed_batches, fh)
        sys.exit(1)


def generate_spectral_lib(config_path: Union[str, Path, Config]):
    """
    Create a SpectralLibrary object and generate the spectral library.

    # TODO full description
    :param config_path: path to config file
    """
    if isinstance(config_path, Config):
        config = config_path
    else:
        config = Config()
        config.read(config_path)
    config.check()

    spec_library = _speclib_from_digestion(config)

    speclib_written_step = ProcessStep(config.output, "speclib_written")
    if not speclib_written_step.is_done():
        results_path = config.output / "results"
        results_path.mkdir(exist_ok=True)

        batch_size = config.speclib_generation_batch_size
        failed_batch_file = config.output / "data" / "speclib_failed_batches.pkl"
        writer, out_file = _get_writer_and_output(results_path, config.output_format)
        batches, mode = _get_batches_and_mode(
            out_file, failed_batch_file, spec_library.obs, batch_size, config.models.get("intensity", "")
        )
        speclib = writer(out_file, mode=mode, min_intensity_threshold=config.min_intensity)
        n_batches = len(batches)

        worker_count = max(1, min(config.num_threads, n_batches))

        with Manager() as manager:
            # setup
            shared_queue = manager.Queue(maxsize=worker_count)
            prediction_progress = manager.Value("i", 0)
            prediction_failure_progress = manager.Value("i", 0)
            writing_progress = manager.Value("i", 0)

            lock = manager.Lock()
            lock_failure = manager.Lock()

            # Create a pool for producer processes
            predictor_pool = pool.Pool(worker_count)

            consumer_process = Process(
                target=speclib.async_write,
                args=(shared_queue, writing_progress, config.custom_to_unimod()),
            )

            try:
                results = []
                for batch in batches:
                    result = predictor_pool.apply_async(
                        _make_predictions,
                        (
                            config,
                            shared_queue,
                            prediction_progress,
                            lock,
                            spec_library.obs.loc[batch],
                        ),
                        error_callback=partial(
                            _make_predictions_error_callback, prediction_failure_progress, lock_failure
                        ),
                    )
                    results.append(result)
                predictor_pool.close()

                with tqdm(
                    total=n_batches, desc="Writing library", postfix={"successful": 0, "missing": 0}
                ) as writer_pbar:
                    consumer_process.start()
                    with tqdm(
                        total=n_batches, desc="Getting predictions", postfix={"successful": 0, "failed": 0}
                    ) as predictor_pbar:
                        while predictor_pbar.n < n_batches:
                            time.sleep(1)
                            pr_fail_val = prediction_failure_progress.value
                            _update(predictor_pbar, {"failed": pr_fail_val, "successful": prediction_progress.value})
                            _update(writer_pbar, {"successful": writing_progress.value, "missing": pr_fail_val})
                        shared_queue.put(None)  # signal the writer process, that it is done
                        predictor_pool.join()  # properly await and terminate the pool

                    while writer_pbar.n < n_batches:  # have to keep updating the writer pbar
                        time.sleep(1)
                        _update(
                            writer_pbar,
                            {"successful": writing_progress.value, "missing": prediction_failure_progress.value},
                        )
                    consumer_process.join()  # properly await the termination of the writer process

                _check_write_failed_batch_file(failed_batch_file, prediction_failure_progress.value, results)

            finally:
                predictor_pool.terminate()
                predictor_pool.join()
                consumer_process.terminate()
                consumer_process.join()

        logger.info("Finished writing the library to disk")
        failed_batch_file.unlink(missing_ok=True)
        speclib_written_step.mark_done()


def _candidate_index(df: pd.DataFrame) -> pd.MultiIndex:
    """Candidate PSMs (raw file, scan number, modified sequence, precursor charge) of a search result / library."""
    return pd.MultiIndex.from_frame(df[pp.SEARCH_RESULT_KEY].astype(str))


def _previous_library(run_dir: Path, spectra_file: Path) -> Optional[Spectra]:
    """
    Annotated and predicted spectra of a spectra file from a previous run, with the original collision energy.

    The collision energy the intensities were predicted at is kept in REUSED_CE_COLUMN; COLLISION_ENERGY is set back
    to the original one of the spectrum (inverse of the run's CE calibration), so the CE calibration of this run applies.

    :param run_dir: output folder of the previous run
    :param spectra_file: the spectra file
    :return: the library, or None if the run lacks the library or CE calibration of the file
    """
    path = run_dir / "data" / spectra_file.with_suffix(".mzml.pred.hdf5").name
    ce_file = run_dir / "results" / "ce_calibration" / f"{spectra_file.stem}_ce.txt"
    if not path.is_file() or not ce_file.is_file():
        logger.warning(f"No library or CE calibration of {spectra_file.stem} in {run_dir}, annotating it again")
        return None
    with open(ce_file) as f:
        best_ce = json.load(f)
    original_ce = {float(ce): float(orig) for orig, ce in best_ce.items()}
    if len(original_ce) != len(best_ce):
        logger.warning(f"CE calibration of {spectra_file.stem} in {run_dir} cannot be inverted, annotating it again")
        return None
    library = Spectra.from_hdf5(path)
    obs = library.obs
    obs[REUSED_CE_COLUMN] = obs["COLLISION_ENERGY"].astype(float)
    is_hcd = (obs["FRAGMENTATION"] == "HCD").to_numpy(dtype=bool)
    original = obs.loc[is_hcd, REUSED_CE_COLUMN].map(original_ce)
    if original.isna().any():
        logger.warning(f"CE of {spectra_file.stem} in {run_dir} does not match its calibration, annotating it again")
        return None
    obs.loc[is_hcd, "COLLISION_ENERGY"] = original
    return library


def _library_from_previous_runs(spectra_file: Path, config: Config) -> Optional[Spectra]:
    """
    Annotated and predicted spectra of the candidates of a spectra file, from previous runs (reuseLibraries).

    Previous runs on the same spectra with the same annotation settings and models, e.g. one run per search engine
    before rescoring the merged search results. Every candidate PSM of this run is taken from the first previous run
    that has it: annotation, predicted iRT and the predicted intensities (predicted again in _reuse_predictions if this
    run's calibrated collision energy differs).

    :param spectra_file: the spectra file
    :param config: the configuration object
    :return: the library, or None (annotate as usual) if a previous run lacks the file or a candidate is not found
    """
    wanted = _candidate_index(pp.load_search(config.output / "msms" / spectra_file.with_suffix(".rescore").name))
    libraries = [_previous_library(run_dir, spectra_file) for run_dir in config.reuse_libraries]
    if any(library is None for library in libraries):
        return None
    if any(not library.var_names.equals(libraries[0].var_names) for library in libraries[1:]):
        logger.warning(f"Previous runs of {spectra_file.stem} annotate different fragment ions, annotating it again")
        return None
    columns = [column for column in libraries[0].obs.columns if all(column in lib.obs for lib in libraries[1:])]
    parts, taken = [], pd.MultiIndex.from_tuples([], names=pp.SEARCH_RESULT_KEY)
    for library in libraries:
        keys = _candidate_index(library.obs)
        keep = np.flatnonzero(keys.isin(wanted) & ~keys.isin(taken) & ~keys.duplicated())
        part = library[keep].copy()
        part.obs = part.obs[columns]
        part.obs_names = [f"{len(parts)}_{i}" for i in range(part.n_obs)]
        parts.append(part)
        taken = taken.append(keys[keep])
    missing = len(wanted.unique()) - len(taken)
    if missing:
        logger.warning(f"{missing} candidates of {spectra_file.stem} not in the previous runs, annotating it again")
        return None
    combined = Spectra(anndata.concat(parts, join="inner", uns_merge="first"))
    combined.obs_names = [str(i) for i in range(combined.n_obs)]
    logger.info(
        f"Reusing the annotation and predictions of {combined.n_obs} candidates of {spectra_file.stem} "
        f"({', '.join(str(len(part)) for part in parts)} from {', '.join(map(str, config.reuse_libraries))})"
    )
    return combined


def _reuse_predictions(library: Spectra, spectra_file: Path, config: Config):
    """
    Keep the reused predictions of a library from previous runs, predicting intensities again where the CE differs.

    Predicted iRT does not depend on the collision energy and is kept for every candidate.

    :param library: library of _library_from_previous_runs after the CE calibration of this run
    :param spectra_file: the spectra file
    :param config: the configuration object
    """
    redo = library.obs[REUSED_CE_COLUMN].to_numpy() != library.obs["COLLISION_ENERGY"].astype(float).to_numpy()
    logger.info(
        f"Reusing the predictions of {int((~redo).sum())} candidates of {spectra_file.stem}, predicting the intensities "
        f"of {int(redo.sum())} at another collision energy"
    )
    if redo.any():
        subset = Spectra(library[np.flatnonzero(redo)].copy())
        subset.obs_names = [str(i) for i in range(subset.n_obs)]
        chunk_idx = None
        if "alphapept" in config.models.get("intensity", "").lower():
            chunk_idx = list(group_iterator(df=subset.obs, group_by_column="PEPTIDE_LENGTH"))
        intensity_predictor = pr.Predictor.from_config(config, model_type="intensity")
        intensity_predictor.predict_intensities(
            data=subset, chunk_idx=chunk_idx, dataset_name=spectra_file.stem, keep_dataset=False
        )
        layer = Spectra.INTENSITY_PRED_LAYER_NAME
        kept, new = np.flatnonzero(~redo), np.flatnonzero(redo)
        stacked = scipy.sparse.vstack([library.layers[layer][kept], subset.layers[layer]]).tocsr()
        library.layers[layer] = stacked[np.argsort(np.concatenate([kept, new]))]
    library.obs.drop(columns=REUSED_CE_COLUMN, inplace=True)


def _ce_calib(spectra_file: Path, config: Config) -> Spectra:
    ce_calib_step = ProcessStep(config.output, "ce_calib." + spectra_file.stem)
    if ce_calib_step.is_done():
        hdf5_path = config.output / "data" / spectra_file.with_suffix(".mzml.pred.hdf5").name
        if hdf5_path.is_file():
            library = Spectra.from_hdf5(hdf5_path)
            return library
        else:
            raise FileNotFoundError(f"{hdf5_path} not found but ce_calib.{spectra_file.stem} found. Please check.")
    tims_meta_file = None
    if config.spectra_type.lower() in ["hdf", "d"]:  # if it is timstof
        tims_meta_file = config.output / "msms" / spectra_file.with_suffix(".timsmeta").name
    aspec = None
    if config.reuse_libraries and not (config.output / "data" / spectra_file.with_suffix(".mzml.hdf5").name).is_file():
        (config.output / "data").mkdir(exist_ok=True)
        aspec = _library_from_previous_runs(spectra_file, config)
    if aspec is None:
        aspec = _annotate_and_get_library(spectra_file, config, tims_meta_file=tims_meta_file)
    _get_best_ce(aspec, spectra_file, config)

    aspec.write_as_hdf5(config.output / "data" / spectra_file.with_suffix(".mzml.pred.hdf5").name)

    ce_calib_step.mark_done()

    return aspec


def run_ce_calibration(
    config_path: Union[str, Path, Config],
):
    """
    Create a CeCalibration object and run the CE calibration.

    # TODO full description
    :param config_path: path to config file
    """
    if isinstance(config_path, Config):
        config = config_path
    else:
        config = Config()
        config.read(config_path)
    config.check()

    # load spectra file names
    spectra_files = pp.list_spectra(input_dir=config.spectra, input_format=config.spectra_type)

    proc_dir = config.output / "proc"
    proc_dir.mkdir(parents=True, exist_ok=True)

    spectra_files = _preprocess(spectra_files, config)

    worker_count = max(1, min(config.num_threads, len(spectra_files)))

    if worker_count > 1:
        processing_pool = JobPool(processes=worker_count)
        for spectra_file in spectra_files:
            processing_pool.apply_async(_ce_calib, [spectra_file, config])
        processing_pool.check_pool()
    else:
        for spectra_file in spectra_files:
            _ce_calib(spectra_file, config)


def _calculate_features(spectra_file: Path, config: Config, xl: bool = False, cms2: bool = False):
    library = _ce_calib(spectra_file, config)
    calc_feature_step_name = "calculate_features_original." if config._run_original else "calculate_features."
    calc_feature_step = ProcessStep(config.output, calc_feature_step_name + spectra_file.stem)
    if calc_feature_step.is_done():
        return

    predict_step = ProcessStep(config.output, "predict." + spectra_file.stem)
    if not config._run_original and not predict_step.is_done():
        if "alphapept" in config.models.get("intensity", "").lower():
            chunk_idx = list(group_iterator(df=library.obs, group_by_column="PEPTIDE_LENGTH"))
        else:
            chunk_idx = None

        intensity_predictor = pr.Predictor.from_config(config, model_type="intensity")
        if xl:
            intensity_predictor.predict_intensities(
                data=library, xl=True, chunk_idx=chunk_idx, dataset_name=spectra_file.stem, keep_dataset=False
            )

            library.write_as_hdf5(config.output / "data" / spectra_file.with_suffix(".mzml.pred.hdf5").name)
            predict_step.mark_done()

        elif REUSED_CE_COLUMN in library.obs:  # library of previous runs (reuseLibraries)
            _reuse_predictions(library, spectra_file, config)
            library.write_as_hdf5(config.output / "data" / spectra_file.with_suffix(".mzml.pred.hdf5").name)
            predict_step.mark_done()

        else:
            intensity_predictor.predict_intensities(
                data=library, chunk_idx=chunk_idx, dataset_name=spectra_file.stem, keep_dataset=False
            )

            irt_predictor = pr.Predictor.from_config(config, model_type="irt")
            irt_predictor.predict_rt(data=library)

            library.write_as_hdf5(config.output / "data" / spectra_file.with_suffix(".mzml.pred.hdf5").name)
            predict_step.mark_done()

    # produce percolator tab files
    fdr_dir = config.output / "results" / config.fdr_estimation_method
    fdr_dir.mkdir(exist_ok=True)
    add_neutral_loss_features = config.ptm_use_neutral_loss
    remove_miss_cleavage_features = ("R" in config.ptm_possible_sites) or ("K" in config.ptm_possible_sites)
    re.generate_features(
        library=library,
        search_type="original",
        output_file=fdr_dir / spectra_file.with_suffix(".original.tab").name,
        additional_columns=config.use_feature_cols,
        all_features=config.all_features,
        xl=xl,
        cms2=cms2,
        regression_method=config.curve_fitting_method,
        featured_ions=None,
    )

    if not config._run_original:
        re.generate_features(
            library=library,
            search_type="rescore",
            output_file=fdr_dir / spectra_file.with_suffix(".rescore.tab").name,
            additional_columns=config.use_feature_cols,
            all_features=config.all_features,
            xl=xl,
            cms2=cms2,
            regression_method=config.curve_fitting_method,
            add_neutral_loss_features=add_neutral_loss_features,
            remove_miss_cleavage_features=remove_miss_cleavage_features,
            task="multifrag" if config.check_multifrag() else "default",
            featured_ions=config.ion_types,
            subtract_peaks=config.chimeric_peak_subtraction,
            search_engine_features=config.search_engine_features,
        )

    calc_feature_step.mark_done()


def _percolator_run(fdr_dir: Path, config: Config, stage: str, num_threads: int, step: ProcessStep, xl: bool = False):
    """
    Run percolator on the original or the rescore features and mark the step as done.

    The rescore features are localized instead if PTMs are localized. The training subset (percolatorSubsetMaxTrain)
    is used for the rescoring runs only, never for the localization models.

    :param fdr_dir: the output directory
    :param config: the configuration object
    :param stage: "original" or "rescore"
    :param num_threads: number of threads used by percolator
    :param step: the process step to mark as done
    :param xl: crosslinked or linear peptide
    """
    features = fdr_dir / f"{stage}.tab"
    if stage == "rescore" and config.ptm_localization:
        re.localize_ptms(
            fdr_dir,
            config.ptm_sites,
            config.ptm_decoy_residues,
            config.ptm_flr_neighbor_split,
            config.ptm_localization_probability,
            config.ptm_decoy_weighting,
            config.ptm_flr_groups,
            config.ptm_flr_estimator,
            config.ptm_sdi_features,
            config.ptm_ambiguity_check,
            config.output / "data",
            num_threads=num_threads,
        )
    elif config.ptm_fdr_control and not xl:
        re.percolator_fdr_modified(
            features, num_threads=num_threads, subset_max_train=config.percolator_subset_max_train
        )
    else:
        re.rescore_with_percolator(
            input_file=features,
            output_folder=fdr_dir,
            xl=xl,
            num_threads=num_threads,
            subset_max_train=config.percolator_subset_max_train,
        )
    if xl:
        output_csms = xl_psm_to_csm(str(fdr_dir), stage, "percolator")
        output_csms = xl_between_or_self(output_csms, score="score")
        xl_preprocessing_plot_csm(str(fdr_dir), output_csms, stage, "percolator")
    step.mark_done()


def _rescore(fdr_dir: Path, config: Config, xl: bool = False):
    """
    High level rescore function for original and rescore.

    :param fdr_dir: the output directory
    :param config: the configuration object
    :param xl: crosslinked or linear peptide
    :raises ValueError: if the provided fdr estimation method in the config is not recognized
    """
    rescore_original_step = ProcessStep(config.output, f"{config.fdr_estimation_method}_original")
    rescore_prosit_step = ProcessStep(config.output, f"{config.fdr_estimation_method}_prosit")

    if config.fdr_estimation_method == "percolator":
        runs = []
        if not rescore_original_step.is_done():
            runs.append(
                lambda threads: _percolator_run(fdr_dir, config, "original", threads, rescore_original_step, xl)
            )
        if not config._run_original and not rescore_prosit_step.is_done():
            runs.append(lambda threads: _percolator_run(fdr_dir, config, "rescore", threads, rescore_prosit_step, xl))
        if len(runs) == 2 and config.num_threads >= 2 and not xl:
            # original and rescore are independent: run them at the same time, sharing numThreads (not for
            # crosslinks, whose csm plots are not thread safe)
            threads = [config.num_threads // 2, config.num_threads - config.num_threads // 2]
            with ThreadPoolExecutor(max_workers=2) as executor:
                for future in [executor.submit(run, n) for run, n in zip(runs, threads, strict=True)]:
                    future.result()
        else:
            for run in runs:
                run(config.num_threads)
    elif config.fdr_estimation_method == "mokapot":
        if not rescore_original_step.is_done():
            re.rescore_with_mokapot(input_file=fdr_dir / "original.tab", output_folder=fdr_dir, xl=xl)
            if xl:
                output_csms_original = xl_psm_to_csm(str(fdr_dir), "original", "mokapot")
                output_csms_original = xl_between_or_self(output_csms_original, score="score")
                xl_preprocessing_plot_csm(str(fdr_dir), output_csms_original, "original", "mokapot")

            rescore_original_step.mark_done()
        if not config._run_original and not rescore_prosit_step.is_done():
            re.rescore_with_mokapot(input_file=fdr_dir / "rescore.tab", output_folder=fdr_dir, xl=xl)
            if xl:
                output_csms_rescore = xl_psm_to_csm(str(fdr_dir), "rescore", "mokapot")
                output_csms_rescore = xl_between_or_self(output_csms_rescore, score="score")
                xl_preprocessing_plot_csm(str(fdr_dir), output_csms_rescore, "rescore", "mokapot")
            rescore_prosit_step.mark_done()
    else:
        raise ValueError(
            'f{config.fdr_estimation_method} is not a valid rescoring tool, use either "percolator" or "mokapot"'
        )


def xl_fdr(df: pd.DataFrame, score: str) -> pd.DataFrame:
    """
    "calculate and add fdr_xl to the DataFrame : (TD-DD)/(TT)".

    :param df: DataFrame containing the data.
    :param score: Column name containing the scores used for calculating FDR
    :return: DataFrame with the new column 'fdr'
    """
    df = df.sort_values(by=score, ascending=False)
    df["TD_sum"] = (df["label"] == "TD").cumsum()
    df["DD_sum"] = (df["label"] == "DD").cumsum()
    df["TT_sum"] = (df["label"] == "TT").cumsum()
    df["q-value"] = (df["TD_sum"] - df["DD_sum"]) / df["TT_sum"]
    df.loc[df["TT_sum"] == 0, "q-value"] = 0  # Handling division by zero
    df = df.drop(["TD_sum", "DD_sum", "TT_sum"], axis=1)
    return df


def xl_between_or_self(df: pd.DataFrame, score: str) -> pd.DataFrame:
    """
    "check csms are self or between links for fdr calculation".

    :param df: DataFrame containing the data.
    :param score: Column name containing the scores used for calculating FDR
    :return: DataFrame with the new column 'fdr'
    """
    df_csms_between = df[df["fdr_group"] == "between"]
    df_csms_self = df[df["fdr_group"] == "self"]
    df_csms_between = xl_fdr(df_csms_between, score=score)
    df_csms_self = xl_fdr(df_csms_self, score=score)
    df_csms = pd.concat([df_csms_between, df_csms_self], axis=0)
    return df_csms


def xl_preprocessing_plot_csm(
    features_dir: str, df: pd.DataFrame, original_or_rescore: str, percolator_or_mokapot: str
):
    """
     Helper function to convert psm to csm level for plotting.

    :param features_dir: the features directory
    :param df: DataFrame containing the data
    :param original_or_rescore: rescoring by search engine features or intensity-based features
    :param percolator_or_mokapot: percolator or mkkapot used for rescoring
    """
    columns_to_keep = [
        "PSMId_a",
        "label",
        "scan_number_a",
        "filename_a",
        "peptide_a",
        "score",
        "q-value",
        "q-value_a",
        "q-value_b",
        "score_a",
        "score_b",
        "fdr_group_a",
        "proteinIds_b",
    ]
    df = df[columns_to_keep]
    df = df.rename(
        columns={
            "PSMId_a": "SpecId",
            "scan_number_a": "ScanNr",
            "filename_a": "filename",
            "peptide_a": "Peptide",
            "fdr_group_a": "fdr_group",
            "proteinIds_b": "Proteins",
        }
    )
    df["is_target"] = df["label"]
    df["label"] = df["label"].replace({"TT": True, "TD": False, "DD": False})
    df_target = df[df["label"]]
    df_decoy = df[~df["label"]]

    if original_or_rescore == "original":
        if percolator_or_mokapot == "percolator":
            df_target.to_csv(features_dir + "/original.percolator.csms.txt", sep="\t", index=False)
            df_decoy.to_csv(features_dir + "/original.percolator.decoy.csms.txt", sep="\t", index=False)
        else:
            df_target.to_csv(features_dir + "/original.mokapot.csms.txt", sep="\t", index=False)
            df_decoy.to_csv(features_dir + "/original.mokapot.decoy.csms.txt", sep="\t", index=False)
    else:
        if percolator_or_mokapot == "percolator":
            df_target.to_csv(features_dir + "/rescore.percolator.csms.txt", sep="\t", index=False)
            df_decoy.to_csv(features_dir + "/rescore.percolator.decoy.csms.txt", sep="\t", index=False)
        else:
            df_target.to_csv(features_dir + "/rescore.mokapot.csms.txt", sep="\t", index=False)
            df_decoy.to_csv(features_dir + "/rescore.mokapot.decoy.csms.txt", sep="\t", index=False)


def xl_psm_to_csm(features_dir: str, original_or_rescore: str, percolator_or_mokapot: str):
    """
     Helper function to convert psm to csm level for ploting.

    :param features_dir: the features directory
    :param original_or_rescore: rescoring by search engine features or intesity-based features
    :param percolator_or_mokapot: percolator or mkkapot used for rescoring
    :return: DataFrame in csm level
    """

    def get_label(row):
        if row["is_decoy_p1_a"] == "False" and row["is_decoy_p2_a"] == "False":
            return "TT"
        elif row["is_decoy_p1_a"] == "True" and row["is_decoy_p2_a"] == "True":
            return "DD"
        else:
            return "TD"

    def min_score(row):
        return min(row["score_a"], row["score_b"])

    if original_or_rescore == "original":
        if percolator_or_mokapot == "percolator":
            decoy_psms = pd.read_csv(features_dir + "/original.percolator.decoy.psms.txt", delimiter="\t")
            target_psms = pd.read_csv(features_dir + "/original.percolator.psms.txt", delimiter="\t")
            psm_id = "PSMId"
        else:
            decoy_psms = pd.read_csv(features_dir + "/original.mokapot.decoy.psms.txt", delimiter="\t")
            target_psms = pd.read_csv(features_dir + "/original.mokapot.psms.txt", delimiter="\t")
            psm_id = "SpecId"

    else:
        if percolator_or_mokapot == "percolator":
            decoy_psms = pd.read_csv(features_dir + "/rescore.percolator.decoy.psms.txt", delimiter="\t")
            target_psms = pd.read_csv(features_dir + "/rescore.percolator.psms.txt", delimiter="\t")
            psm_id = "PSMId"
        else:
            decoy_psms = pd.read_csv(features_dir + "/rescore.mokapot.decoy.psms.txt", delimiter="\t")
            target_psms = pd.read_csv(features_dir + "/rescore.mokapot.psms.txt", delimiter="\t")
            psm_id = "SpecId"

    new_columns = [
        "raw_file",
        "scan_number",
        "mod_pep_a",
        "mod_pep_b",
        "charge",
        "decoy_p1",
        "is_decoy_p1",
        "decoy_p2",
        "is_decoy_p2",
        "fdr_group",
        "base_sequence_p1",
        "base_sequence_p2",
        "index",
    ]
    split_data = target_psms[psm_id].str.rsplit("-", n=12, expand=True)
    split_data.columns = new_columns
    df_psm_target = pd.concat([target_psms, split_data], axis=1)
    split_data = decoy_psms[psm_id].str.rsplit("-", n=12, expand=True)
    split_data.columns = new_columns
    df_psm_decoy = pd.concat([decoy_psms, split_data], axis=1)
    df_psm = pd.concat([df_psm_decoy, df_psm_target], axis=0)
    df_psm[
        [
            "index_csm",
            "_",
            "which_pep",
        ]
    ] = df_psm["index"].str.split("_", expand=True)
    df_psm.drop(columns=["index", "_", "decoy_p1", "decoy_p2"], inplace=True)
    df_pep_1 = df_psm[df_psm["which_pep"] == "1"].copy()
    df_pep_2 = df_psm[df_psm["which_pep"] == "2"].copy()
    df_pep_1.drop(columns=["which_pep"], inplace=True)
    df_pep_2.drop(columns=["which_pep"], inplace=True)
    df_pep_1.columns = [col + "_a" if col != "index_csm" else col for col in df_pep_1.columns]
    df_pep_2.columns = [col + "_b" if col != "index_csm" else col for col in df_pep_2.columns]
    df_csm = pd.merge(df_pep_1, df_pep_2, on="index_csm")
    df_csm["score"] = df_csm.apply(min_score, axis=1)
    df_csm["label"] = df_csm.apply(get_label, axis=1)
    df_csm.rename(columns={"fdr_group_b": "fdr_group"}, inplace=True)
    return df_csm


def prepare_rescore_xl_psm_level(features_dir: str, original_or_rescore: str):
    """
    Helper function for running percolator on psm level instead of csm level.

    :param features_dir: the features directory
    :param original_or_rescore: rescoring by search engine features or intensity-based features
    :return: DataFrame on psm level
    """

    def extract_label_pep_a(specid):
        if "-decoy_p1-True-" in specid:
            return -1
        elif "-decoy_p1-False-" in specid:
            return 1
        else:
            return None

    def extract_label_pep_b(specid):
        if "-decoy_p2-True-" in specid:
            return -1
        elif "-decoy_p2-False-" in specid:
            return 1
        else:
            return None

    columns_to_remove_psm_rescore = [
        "run_name",
        "scan_number",
        "precursor_mass",
        "precursor_charge",
        "crosslinker_name",
        "aa_len_p1",
        "link_pos_p1",
        "linked_aa_p1",
        "mods_p1",
        "sequence_p1",
        "sequence_p2",
        "start_pos_p1",
        "start_pos_p2",
        "mod_pos_p1",
        "aa_len_p2",
        "link_pos_p2",
        "linked_aa_p2",
        "mods_p2",
        "mod_pos_p2",
        "linear",
        "match_score",
        "decoy",
        "RAW_FILE",
        "MASS",
        "PRECURSOR_CHARGE",
        "CROSSLINKER_TYPE",
        "REVERSE",
        "SCAN_NUMBER",
        "SEQUENCE_A",
        "SEQUENCE_B",
        "Modifications_A",
        "Modifications_B",
        "CROSSLINKER_POSITION_A",
        "CROSSLINKER_POSITION_B",
        "ModificationPositions1",
        "MZ_RANGE",
        "ModificationPositions2",
        "MODIFIED_SEQUENCE_A",
        "MODIFIED_SEQUENCE_B",
        "FRAGMENTATION",
        "INSTRUMENT_TYPES",
        "MASS_ANALYZER",
        "Proteins",
        "spectral_angle",
        "SCORE",
        "Peptide Position 1",
        "Peptide Position 2",
        "Protein 1",
        "Protein 2",
        "Unnamed: 0",
    ]
    columns_to_remove_psm_original = list(set(columns_to_remove_psm_rescore + ["match_score"]) - {"spectral_angle"})

    if original_or_rescore == "original":
        columns_to_remove_psm = columns_to_remove_psm_original
        rescore_tab_file = pd.read_csv(features_dir + "/original.tab", sep="\t")
    else:
        columns_to_remove_psm = columns_to_remove_psm_rescore
        rescore_tab_file = pd.read_csv(features_dir + "/rescore.tab", sep="\t")

    rescore_tab_file.drop(columns=columns_to_remove_psm, inplace=True, errors="ignore")
    rescore_tab_file = rescore_tab_file.fillna(0)
    string_columns = [
        "SpecId",
        "Peptide",
        "Proteins",
        "protein_p1",
        "protein_p2",
        "filename",
        "fdr_group",
        "base_sequence_p1",
        "base_sequence_p2",
    ]
    rescore_tab_file_numeric_columns = [col for col in rescore_tab_file.columns if col not in string_columns]
    rescore_tab_file[rescore_tab_file_numeric_columns] = rescore_tab_file[rescore_tab_file_numeric_columns].apply(
        pd.to_numeric, errors="coerce"
    )
    rescore_tab_file = rescore_tab_file.reset_index(drop=True)
    rescore_tab_file["Proteins"] = (
        "p1_" + rescore_tab_file["protein_p1"].astype(str) + "_p2_" + rescore_tab_file["protein_p2"].astype(str)
    )
    rescore_tab_file["SpecId"] = (
        rescore_tab_file["SpecId"]
        + "-decoy_p1-"
        + rescore_tab_file["decoy_p1"].astype(str)
        + "-decoy_p2-"
        + rescore_tab_file["decoy_p2"].astype(str)
        + "-"
        + rescore_tab_file["fdr_group"].astype(str)
        + "-"
        + rescore_tab_file["base_sequence_p1"].astype(str)
        + "-"
        + rescore_tab_file["base_sequence_p2"].astype(str)
        + "-"
        + rescore_tab_file.index.astype(str)  # Adding index to SpecId
    )
    rescore_tab_file.drop(
        columns=[
            "protein_p1",
            "protein_p2",
            "decoy_p1",
            "decoy_p2",
            "fdr_group",
            "base_sequence_p1",
            "base_sequence_p2",
        ],
        inplace=True,
    )
    columns_feature_pep_a = [
        col
        for col in rescore_tab_file.columns
        if col.endswith("_a") or col.endswith("_A") or not (col.endswith("_b") or col.endswith("_B"))
    ]
    columns_feature_pep_b = [
        col
        for col in rescore_tab_file.columns
        if col.endswith("_b") or col.endswith("_B") or not (col.endswith("_a") or col.endswith("_A"))
    ]
    rescore_tab_file_a = rescore_tab_file.loc[:, columns_feature_pep_a]
    rescore_tab_file_b = rescore_tab_file.loc[:, columns_feature_pep_b]
    rescore_tab_file_a["SpecId"] = rescore_tab_file_a["SpecId"] + "_" + rescore_tab_file_a.index.astype(str)
    rescore_tab_file_b["SpecId"] = rescore_tab_file_b["SpecId"] + "_" + rescore_tab_file_b.index.astype(str)
    # 1 means pep a and 2 means pep b
    rescore_tab_file_a["SpecId"] = rescore_tab_file_a["SpecId"] + "_1"
    rescore_tab_file_b["SpecId"] = rescore_tab_file_b["SpecId"] + "_2"
    rescore_tab_file_a.columns = [col[:-2] if col.endswith(("_a", "_A")) else col for col in rescore_tab_file_a.columns]
    rescore_tab_file_b.columns = [col[:-2] if col.endswith(("_b", "_B")) else col for col in rescore_tab_file_b.columns]
    rescore_tab_file_a["Label"] = rescore_tab_file_a["SpecId"].apply(extract_label_pep_a)
    rescore_tab_file_b["Label"] = rescore_tab_file_b["SpecId"].apply(extract_label_pep_b)
    string_columns = ["SpecId", "Peptide", "Proteins", "filename"]
    rescore_tab_file_numeric_columns = [col for col in rescore_tab_file_a.columns if col not in string_columns]
    rescore_tab_file_a[rescore_tab_file_numeric_columns] = rescore_tab_file_a[rescore_tab_file_numeric_columns].apply(
        pd.to_numeric, errors="coerce"
    )
    rescore_tab_file_b[rescore_tab_file_numeric_columns] = rescore_tab_file_b[rescore_tab_file_numeric_columns].apply(
        pd.to_numeric, errors="coerce"
    )
    # change ExpMass of rescore_tab_file_a
    max_expmass = rescore_tab_file_a["ExpMass"].max()
    rescore_tab_file_b["ExpMass"] += max_expmass
    input_psm_rescore = pd.concat([rescore_tab_file_a, rescore_tab_file_b], axis=0, ignore_index=True)
    input_psm_rescore["Proteins"].fillna("unknown", inplace=True)
    return input_psm_rescore


def input_xifdr(fdr_dir: str, xisearch_or_scout: str):
    """
     Helper function to generate xiFDR input.

    :param fdr_dir: the output directory
    :param xisearch_or_scout: applied search engine
    """

    def convert_percolator_output(df: pd.DataFrame):
        df["SpecId_raw_name_scan"] = df["SpecId"].str.extract(r"^([^-]+-[^-]+)")
        split_data = df["SpecId"].str.rsplit("-", n=12, expand=True)
        new_columns = [
            "raw_file",
            "scan_number",
            "mod_pep_a",
            "mod_pep_b",
            "charge",
            "decoy_p1",
            "is_decoy_p1",
            "decoy_p2",
            "is_decoy_p2",
            "fdr_group",
            "base_sequence_p1",
            "base_sequence_p2",
            "index",
        ]
        split_data.columns = new_columns
        df = pd.concat([df, split_data], axis=1)
        df.drop(columns=["index", "decoy_p1", "decoy_p2"], inplace=True)
        df["SpecId"] = df["SpecId"].str.split("-decoy_p1").str[0]
        df = df.reset_index()
        df["protein_p1"] = df["Proteins"].apply(lambda x: x.split("_p2_")[0].split("p1_")[1])
        df["protein_p2"] = df["Proteins"].apply(lambda x: x.split("_p2_")[1])
        return df

    def input_columns_xifdr(df: pd.DataFrame, input_type: str):
        new_column_names = {
            "run_name": "run",
            "FileName": "run",
            "scan_number": "scan",
            "ScanNumber": "scan",
            "sequence_p1": "peptide1",
            "MODIFIED_SEQUENCE_A_xiFDR": "peptide1",
            "sequence_p2": "peptide2",
            "MODIFIED_SEQUENCE_B_xiFDR": "peptide2",
            "link_pos_p1": "peptide link 1",
            "AlphaPos": "peptide link 1",
            "link_pos_p2": "peptide link 2",
            "BetaPos": "peptide link 2",
            "is_decoy_p1": "is decoy 1",
            "is_decoy_p2": "is decoy 2",
            "precursor_charge": "precursor charge",
            "Charge": "precursor charge",
            "protein_p1": "accession1",
            "Protein 1": "accession1",
            "protein_p2": "accession2",
            "Protein 2": "accession2",
            "start_pos_p1": "peptide position 1",
            "Peptide Position 1": "peptide position 1",
            "start_pos_p2": "peptide position 2",
            "Peptide Position 2": "peptide position 2",
        }
        if input_type == "xisearch":
            new_column_names["match_score"] = "score"
        elif input_type == "percolator":
            new_column_names["score_percolator"] = "score"
        elif input_type == "scout":
            new_column_names["ClassificationScore"] = "score"

        df = df.rename(columns=new_column_names)

        return df

    xifdr_columns = [
        "run",
        "scan",
        "peptide1",
        "peptide2",
        "peptide link 1",
        "peptide link 2",
        "is decoy 1",
        "is decoy 2",
        "precursor charge",
        "accession1",
        "accession2",
        "peptide position 1",
        "peptide position 2",
        "score",
    ]

    # covert percolator output for merging with internal search engine output
    df_percolator_csm_target = pd.read_csv(str(fdr_dir) + "/rescore.percolator.csms.txt", sep="\t")
    df_percolator_csm_target = convert_percolator_output(df_percolator_csm_target)
    df_percolator_csm_decoy = pd.read_csv(str(fdr_dir) + "/rescore.percolator.decoy.csms.txt", sep="\t")
    df_percolator_csm_decoy = convert_percolator_output(df_percolator_csm_decoy)
    df_percolator_csm = pd.concat([df_percolator_csm_target, df_percolator_csm_decoy], ignore_index=True)
    df_percolator_csm.rename(columns={"score": "score_percolator"}, inplace=True)
    df_percolator_csm.reset_index(drop=True, inplace=True)

    # read df_search_engine_internal
    df_search_engine_internal = pd.read_csv(str(fdr_dir) + "/rescore_features_csm.tab", sep="\t")
    df_search_engine_internal.drop("SCORE", axis=1, inplace=True)
    if xisearch_or_scout == "scout":
        df_search_engine_internal["MODIFIED_SEQUENCE_A_xiFDR"] = df_search_engine_internal["AlphaPeptide"].str.replace(
            "+", ""
        )
        df_search_engine_internal["MODIFIED_SEQUENCE_B_xiFDR"] = df_search_engine_internal["BetaPeptide"].str.replace(
            "+", ""
        )
    df_search_engine_internal.reset_index(drop=True, inplace=True)

    # merge percolator and search engine output
    merged_xisearch_percolator = pd.merge(df_search_engine_internal, df_percolator_csm, on="SpecId")
    merged_xisearch_percolator = merged_xisearch_percolator.rename(columns=lambda x: x[:-2] if x.endswith("_x") else x)
    df_percolator_xifdr_input = input_columns_xifdr(merged_xisearch_percolator, input_type="percolator")
    df_percolator_xifdr_input = df_percolator_xifdr_input[xifdr_columns]
    df_percolator_xifdr_input = df_percolator_xifdr_input.loc[:, ~df_percolator_xifdr_input.columns.duplicated()]
    df_percolator_xifdr_input.reset_index(drop=True, inplace=True)
    df_percolator_xifdr_input.to_csv(str(fdr_dir) + "/percolator_xifdr_input.csv", sep=",", index=False)

    if xisearch_or_scout == "xisearch":
        df_xisearch_xifdr_input = input_columns_xifdr(merged_xisearch_percolator, input_type="xisearch")
        df_xisearch_xifdr_input = df_xisearch_xifdr_input[xifdr_columns]
        df_xisearch_xifdr_input = df_xisearch_xifdr_input.loc[:, ~df_xisearch_xifdr_input.columns.duplicated()]
        df_xisearch_xifdr_input.reset_index(drop=True, inplace=True)
        df_xisearch_xifdr_input.to_csv(str(fdr_dir) + "/xisearch_xifdr_input.csv", sep=",", index=False)
    elif xisearch_or_scout == "scout":
        df_scout_xifdr_input = input_columns_xifdr(merged_xisearch_percolator, input_type="scout")
        df_scout_xifdr_input = df_scout_xifdr_input[xifdr_columns]
        df_scout_xifdr_input = df_scout_xifdr_input.loc[:, ~df_scout_xifdr_input.columns.duplicated()]
        df_scout_xifdr_input.reset_index(drop=True, inplace=True)
        df_scout_xifdr_input.to_csv(str(fdr_dir) + "/scout_xifdr_input.csv", sep=",", index=False)


def run_rescoring(config_path: Union[str, Path, Config]):
    """
    Create a ReScore object and run the rescoring.

    # TODO full description
    :param config_path: path to config file
    """
    if isinstance(config_path, Config):
        config = config_path
    else:
        config = Config()
        config.read(config_path)
    config.check()

    # load spectra file names
    spectra_files = pp.list_spectra(input_dir=config.spectra, input_format=config.spectra_type)

    proc_dir = config.output / "proc"
    proc_dir.mkdir(parents=True, exist_ok=True)

    spectra_files = _preprocess(spectra_files, config)

    worker_count = max(1, min(config.num_threads, len(spectra_files)))

    if worker_count > 1:
        processing_pool = JobPool(processes=worker_count)
        for spectra_file in spectra_files:
            _ = processing_pool.apply_async(_ce_calib, [spectra_file, config])
        processing_pool.check_pool()
    else:
        for spectra_file in spectra_files:
            _ = _ce_calib(spectra_file, config)

    if worker_count > 1:
        processing_pool = JobPool(processes=worker_count)
        for spectra_file in spectra_files:
            if "xl" in config.models.get("intensity", "").lower():
                if "cms2" in config.models.get("intensity", "").lower():
                    cms2 = True
                else:
                    cms2 = False
                processing_pool.apply_async(_calculate_features, [spectra_file, config], xl=True, cms2=cms2)
            else:
                processing_pool.apply_async(_calculate_features, [spectra_file, config])
        processing_pool.check_pool()
    else:
        for spectra_file in spectra_files:
            if "xl" in config.models.get("intensity", "").lower():
                if "cms2" in config.models.get("intensity", "").lower():
                    cms2 = True
                else:
                    cms2 = False
                _calculate_features(spectra_file, config, xl=True, cms2=cms2)
            else:
                _calculate_features(spectra_file, config)

    # prepare rescoring
    fdr_dir = config.output / "results" / config.fdr_estimation_method
    original_tab_files = [fdr_dir / spectra_file.with_suffix(".original.tab").name for spectra_file in spectra_files]

    prepare_tab_original_step = ProcessStep(config.output, f"{config.fdr_estimation_method}_prepare_tab_original")

    if not prepare_tab_original_step.is_done():
        logger.info("Merging input tab files for rescoring without peptide property prediction")
        re.merge_input(
            tab_files=original_tab_files, output_file=fdr_dir / "original.tab", num_threads=config.num_threads
        )
        prepare_tab_original_step.mark_done()

    if not config._run_original:
        rescore_tab_files = [fdr_dir / spectra_file.with_suffix(".rescore.tab").name for spectra_file in spectra_files]
        prepare_tab_rescore_step = ProcessStep(config.output, f"{config.fdr_estimation_method}_prepare_tab_prosit")
        if not prepare_tab_rescore_step.is_done():
            logger.info("Merging input tab files for rescoring with peptide property prediction")
            re.merge_input(
                tab_files=rescore_tab_files, output_file=fdr_dir / "rescore.tab", num_threads=config.num_threads
            )
            prepare_tab_rescore_step.mark_done()

    # rescoring (the original features are always rescored when generating training data)
    perform_rescoring = config.perform_rescoring or config._run_original
    if "xl" in config.models.get("intensity", "").lower():  # xl-psm-level
        rescore_features_path = fdr_dir / "rescore_features_csm.tab"
        if not rescore_features_path.exists():
            shutil.copy(fdr_dir / "rescore.tab", rescore_features_path)
            input_psm_rescore = prepare_rescore_xl_psm_level(str(fdr_dir), "rescore")
            input_psm_rescore.to_csv(str(fdr_dir) + "/rescore.tab", sep="\t", index=None)

        original_features_path = fdr_dir / "original_features_csm.tab"
        if not config._run_original and not original_features_path.exists():
            shutil.copy(fdr_dir / "original.tab", original_features_path)
            input_psm_original = prepare_rescore_xl_psm_level(str(fdr_dir), "original")
            input_psm_original.to_csv(str(fdr_dir) + "/original.tab", sep="\t", index=None)

        if not perform_rescoring:
            logger.info("Skipping rescoring, as perform_rescoring is set to false.")
            return

        _rescore(fdr_dir, config, xl=True)

        logger.info("Finished rescoring.")
        if not config._run_original:
            generate_xifdr_input_step = ProcessStep(config.output, "generate_xifdr_input")
            if not generate_xifdr_input_step.is_done():
                logger.info("Generating xiFDR input.")
                if config.search_results_type == "xisearch":
                    input_xifdr(str(fdr_dir), "xisearch")
                elif config.search_results_type == "scout":
                    input_xifdr(str(fdr_dir), "scout")
                logger.info("Finished Generating xiFDR input.")
                generate_xifdr_input_step.mark_done()

    else:
        if not perform_rescoring:
            logger.info("Skipping rescoring, plotting and quantification, as perform_rescoring is set to false.")
            return

        _rescore(fdr_dir, config)
        # plotting
        if not config.perform_plotting:
            logger.info("Skipping summary plots, as perform_plotting is set to false.")
        else:
            logger.info("Generating summary plots...")
            if not config.ptm_localization and not config._run_original:
                pl.plot_all(fdr_dir, config)
                mod_fdr_dir = fdr_dir / "mod_fdr"
                if (
                    config.ptm_fdr_control
                    and config.fdr_estimation_method == "percolator"
                    and (mod_fdr_dir / "rescore.percolator.psms.txt").is_file()
                ):
                    pl.plot_all(mod_fdr_dir, config)
        logger.info("Finished rescoring.")

        if config.quantification:
            logger.info("Starting quantification")
            # method contains picked-group-FDR call
            apply_quant(config)
            logger.info("Finished quantification")


def run_generate_training_data(config_path: Union[str, Path, Config]):
    """Generate training data for Prosit from the rescoring results. This function will run rescoring with the original features and then export the training data as parquet files."""
    if isinstance(config_path, Config):
        config = config_path
    else:
        config = Config()
        config.read(config_path)
    config.check()

    if config.fdr_estimation_method != "percolator":
        raise ValueError("GenerateTrainingData requires fdr_estimation_method to be 'percolator'.")

    try:
        config._run_original = True
        run_rescoring(config)
    finally:
        config._run_original = False

    _export_training_data_parquets(config=config)


def _export_training_data_parquets(config: Config, fdr_threshold: float = 0.01):
    """
    Export full and original Percolator 1% FDR filtered training data parquet files.

    :param config: Parsed Oktoberfest configuration
    :param fdr_threshold: Strict q-value threshold used for the filtered output
    """
    if config.fdr_estimation_method != "percolator":
        raise ValueError("GenerateTrainingData parquet export only supports Percolator output.")

    data_dir = config.output / "data"
    fdr_dir = config.output / "results" / "percolator"
    original_psms_path = fdr_dir / "original.percolator.psms.txt"
    original_tab_path = fdr_dir / "original.tab"

    if not original_psms_path.is_file():
        raise FileNotFoundError(f"Expected Percolator target PSM output not found: {original_psms_path}")
    if not original_tab_path.is_file():
        raise FileNotFoundError(f"Expected Percolator original tab file not found: {original_tab_path}")

    original_psms = pd.read_csv(original_psms_path, sep="\t")
    _require_columns(original_psms, ["PSMId", "q-value"], original_psms_path)
    original_psms["q-value"] = pd.to_numeric(original_psms["q-value"])
    accepted_psm_ids = set(original_psms.loc[original_psms["q-value"] < fdr_threshold, "PSMId"].astype(str))

    original_tab = pd.read_csv(original_tab_path, sep="\t")
    _require_columns(original_tab, ["SpecId", "filename"], original_tab_path)
    original_tab["SpecId"] = original_tab["SpecId"].astype(str)
    original_tab["filename"] = original_tab["filename"].astype(str)

    for run_name in original_tab["filename"].drop_duplicates():
        hdf5_path = data_dir / f"{run_name}.mzml.hdf5"
        if not hdf5_path.is_file():
            raise FileNotFoundError(f"Expected annotated spectra HDF5 not found: {hdf5_path}")

        full_parquet_path = data_dir / f"{run_name}.parquet"
        fdr_parquet_path = data_dir / f"{run_name}.fdr1.parquet"

        logger.info(f"Writing full training parquet file to {full_parquet_path}")
        pp.convert_anndata_to_parquet(hdf5_path, full_parquet_path)

        library = Spectra.from_hdf5(hdf5_path)
        run_tab = original_tab[original_tab["filename"] == run_name]
        if len(run_tab) != library.n_obs:
            raise ValueError(
                f"Per-run original.tab row count for {run_name!r} ({len(run_tab)}) does not match "
                f"HDF5 row count ({library.n_obs}) in {hdf5_path}."
            )

        accepted_run_rows = run_tab["SpecId"].isin(accepted_psm_ids).to_numpy()
        filtered_library = Spectra(library[accepted_run_rows].copy())
        with tempfile.TemporaryDirectory(dir=data_dir) as tmp_dir:
            filtered_hdf5_path = Path(tmp_dir) / f"{run_name}.fdr1.hdf5"
            filtered_library.write_as_hdf5(filtered_hdf5_path)
            logger.info(f"Writing 1% FDR training parquet file to {fdr_parquet_path}")
            pp.convert_anndata_to_parquet(filtered_hdf5_path, fdr_parquet_path)


def _require_columns(df: pd.DataFrame, columns: list[str], file_path: Path):
    missing_columns = [column for column in columns if column not in df.columns]
    if missing_columns:
        raise ValueError(f"{file_path} is missing required column(s): {', '.join(missing_columns)}")


def _extract_gene_name(protein: str) -> Optional[str]:
    """Extract the gene name from a UniProt-style FASTA header.

    Expects the ``db|accession|entry`` format (e.g. ``sp|P16402|H13_HUMAN``) and returns
    the gene portion of the entry name (``H13`` for the example). Returns ``None`` when the
    header does not follow this format or no gene name can be extracted.

    :param protein: single protein identifier / FASTA header
    :return: extracted gene name, or ``None`` if it cannot be determined
    """
    parts = protein.split("|")
    if len(parts) != 3:
        return None
    entry_name = parts[2]
    gene = entry_name.split("_")[0]
    return gene or None


def _build_speclib_rows_for_file(filename: str, peptide_df: pd.DataFrame, data_dir: Path) -> list[list]:
    """Build DIA spectral library rows for a single raw file.

    :param filename: raw file name (without the .mzml.pred.hdf5 suffix)
    :param peptide_df: percolator PSMs belonging to this file
    :param data_dir: directory holding the annotated ``*.mzml.pred.hdf5`` files
    :raises FileNotFoundError: if the annotated spectra HDF5 for the file is missing
    :return: list of unexploded spectral library rows following ``_COLUMNS_DIA_SPECLIB``
    """
    hdf5_path = data_dir / f"{filename}.mzml.pred.hdf5"
    if not hdf5_path.is_file():
        raise FileNotFoundError(f"Expected annotated spectra HDF5 not found: {hdf5_path}")

    library = Spectra.from_hdf5(hdf5_path)
    obs = library.obs

    spec_id_cols = ["RAW_FILE", "SCAN_NUMBER", "MODIFIED_SEQUENCE", "PRECURSOR_CHARGE"]
    if "SCAN_EVENT_NUMBER" in obs.columns:
        spec_id_cols.append("SCAN_EVENT_NUMBER")
    obs["PSMId"] = obs[spec_id_cols].astype(str).agg("-".join, axis=1)

    # Filter 1% FDR peptide level
    obs["avg_retention_time_sec"] = obs.groupby(["MODIFIED_SEQUENCE", "PRECURSOR_CHARGE"])[
        "retention_time_sec"
    ].transform("mean")

    # PSMId only ends with the charge if there is no SCAN_EVENT_NUMBER, so take the charge from the spectra
    peptide_df = peptide_df.copy()
    peptide_df["charge"] = peptide_df["PSMId"].map(obs.set_index("PSMId")["PRECURSOR_CHARGE"])
    peptide_df.drop_duplicates(["peptide", "charge"], inplace=True)
    accepted_peptide_ids = peptide_df.loc[peptide_df["q-value"] <= 0.01, "PSMId"]
    positions = np.flatnonzero(obs["PSMId"].isin(accepted_peptide_ids).to_numpy())

    mz_layer = library.layers["mz"]
    int_layer = library.layers["pred_int"]
    annotations = library.var.index.to_numpy()
    ion_types = library.var["type"].to_numpy()
    frag_charges = library.var["charge"].to_numpy()
    frag_nums = library.var["num"].to_numpy()

    meta = obs.iloc[positions][
        [
            "MASS",
            "PROTEINS",
            "RAW_FILE",
            "SEQUENCE",
            "MODIFIED_SEQUENCE",
            "PRECURSOR_CHARGE",
            "RETENTION_TIME",
            "avg_retention_time_sec",
        ]
    ]

    # precursor m/z, same definition as the ExpMass percolator feature
    precursor_mz = (meta["MASS"] + meta["PRECURSOR_CHARGE"] * c.PARTICLE_MASSES["PROTON"]) / meta["PRECURSOR_CHARGE"]

    rows = []
    for pos, row, mz_precursor in zip(positions, meta.itertuples(index=False), precursor_mz, strict=True):
        mz = mz_layer[pos].toarray().ravel()
        intensity = int_layer[pos].toarray().ravel()
        keep = (mz > 0.0) & (intensity > 0.05)

        proteins = row.PROTEINS
        protein_list = proteins.split(";")
        proteotypic = 1 if len(protein_list) == 1 else 0

        gene_list = [_extract_gene_name(protein) for protein in protein_list]
        genes = ";".join(str(gene) for gene in gene_list)

        rows.append(
            [
                mz_precursor,
                mz[keep],
                annotations[keep].tolist(),
                protein_list[0],
                gene_list[0],
                row.SEQUENCE,
                row.MODIFIED_SEQUENCE,
                row.PRECURSOR_CHARGE,
                intensity[keep] * 10000,
                row.RETENTION_TIME,
                None,  # PrecursorIonMobility
                ion_types[keep].tolist(),
                frag_charges[keep].tolist(),
                frag_nums[keep].tolist(),
                None,  # FragmentLossType
                row.avg_retention_time_sec,
                proteins,
                genes,
                proteotypic,
            ]
        )
    return rows


def generate_spectral_lib_fdr_control(config_path: Union[str, Path, Config]):
    """Generate spectral library from experimental data.

    :param config_path: Path to the configuration file.
    :raises ValueError: If the FDR estimation method is not Percolator.
    :raises FileNotFoundError: If the annotated spectra HDF5 for a file is missing.
    """
    if isinstance(config_path, Config):
        config = config_path
    else:
        config = Config()
        config.read(config_path)
    config.check()

    if config.fdr_estimation_method != "percolator":
        raise ValueError("GenerateTrainingData parquet export only supports Percolator output.")

    data_dir = config.output / "data"
    fdr_dir = config.output / "results" / "percolator"
    peptide_path = fdr_dir / "rescore.percolator.psms.txt"

    peptide_df = pd.read_csv(peptide_path, sep="\t")
    peptide_df["q-value"] = pd.to_numeric(peptide_df["q-value"], errors="coerce")
    grouped_files = [(filename, group) for filename, group in peptide_df.groupby("filename")]

    if config.num_threads > 1:
        processing_pool = JobPool(processes=config.num_threads)
        for filename, group in grouped_files:
            processing_pool.apply_async(_build_speclib_rows_for_file, [filename, group, data_dir])
        rows = [row for file_rows in processing_pool.check_pool() for row in file_rows]
    else:
        rows = []
        for filename, group in grouped_files:
            rows.extend(_build_speclib_rows_for_file(filename, group, data_dir))

    df = pd.DataFrame(
        rows,
        columns=[
            "PrecursorMz",
            "ProductMz",
            "Annotation",
            "ProteinId",
            "GeneName",
            "PeptideSequence",
            "ModifiedPeptideSequence",
            "PrecursorCharge",
            "LibraryIntensity",
            "NormalizedRetentionTime",
            "PrecursorIonMobility",
            "FragmentType",
            "FragmentCharge",
            "FragmentSeriesNumber",
            "FragmentLossType",
            "AverageExperimentalRetentionTime",
            "AllMappedProteins",
            "AllMappedGenes",
            "Proteotypic",
        ],
    )
    fragment_cols = [
        "ProductMz",
        "Annotation",
        "LibraryIntensity",
        "FragmentType",
        "FragmentCharge",
        "FragmentSeriesNumber",
    ]
    df = df.explode(fragment_cols, ignore_index=True)
    # keep the most intense entry per fragment
    df = df.sort_values(by="LibraryIntensity", ascending=False)
    df = df.drop_duplicates(["ModifiedPeptideSequence", "PrecursorCharge", "ProductMz", "NormalizedRetentionTime"])
    df.sort_values(by=["ModifiedPeptideSequence", "PrecursorCharge", "ProductMz"], inplace=True)
    df.to_csv(config.output / "results/spec_lib.tsv", sep="\t", index=False)


def run_job(config_path: Union[str, Path, Config]):
    """
    Run oktoberfest based on job type given in the config file.

    :param config_path: Path to config file
    :raises ValueError: In case the job_type in the provided config file is not known
    """
    if isinstance(config_path, Config):
        config = config_path
    else:
        config = Config()
        config.read(config_path)
    config.check()

    output_folder = config.output
    job_type = config.job_type

    output_folder.mkdir(exist_ok=True)

    # add file handler to root logger
    base_logger = logging.getLogger()
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(name)s::%(funcName)s %(message)s")
    suffix = datetime.datetime.now().strftime("%Y_%m_%d-%H_%M_%S")
    logging_output = output_folder / f"{job_type}_{suffix}.log"
    file_handler = logging.FileHandler(filename=logging_output)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    base_logger.addHandler(file_handler)

    logger.info(f"Oktoberfest version {__version__}\n{__copyright__}")
    logger.info("Job executed with the following config:")
    logger.info(json.dumps(config.data, indent=4))

    try:
        if job_type == "SpectralLibraryGeneration":
            generate_spectral_lib(config_path)
        elif job_type == "GenerateTrainingData":
            run_generate_training_data(config_path)
        elif job_type == "CollisionEnergyCalibration":
            run_ce_calibration(config_path)
        elif job_type == "Rescoring":
            logger.info("Starting rescoring")
            run_rescoring(config_path)
            if config.acquisition_method == "dia" and config.perform_rescoring:
                generate_spectral_lib_fdr_control(config_path)
        else:
            raise ValueError(f"Unknown job_type in config: {job_type}")
    finally:
        file_handler.close()
        base_logger.removeHandler(file_handler)
