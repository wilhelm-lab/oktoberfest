import datetime
import json
import logging
import pickle
import sys
import time
from functools import partial
from math import ceil
from multiprocessing import Manager, Process, pool
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd
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

logger = logging.getLogger(__name__)


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
        if not config.search_results_type == "internal":
            logger.info(f"Converting search results from {config.search_results} to internal search result.")

            msms_output = config.output / "msms"
            msms_output.mkdir(exist_ok=True)
            internal_search_file = msms_output / "msms.prosit"
            tmt_label = config.tag
            ptm_unimods = config.ptm_unimod_id
            ptm_sites = config.ptm_possible_sites
            search_results = pp.convert_search(
                input_path=config.search_results,
                search_engine=config.search_results_type,
                tmt_label=tmt_label,
                custom_mods=config.custom_to_unimod(),
                output_file=internal_search_file,
                ptm_unimod_id=ptm_unimods,
                ptm_sites=ptm_sites,
            )
            if config.spectra_type.lower() in ["d", "hdf"]:
                timstof_metadata = pp.convert_timstof_metadata(
                    input_path=config.search_results,
                    search_engine=config.search_results_type,
                    output_file=msms_output / "tims_meta.csv",
                )
        else:
            internal_search_file = config.search_results
            search_results = pp.load_search(internal_search_file)

            # TODO add support for internal timstof metadata
        logger.info(f"Read {len(search_results)} PSMs from {internal_search_file}")

        # filter search results
        if config.predict_intensity_locally:
            model_type = Path(config.models["intensity"]).stem
        else:
            model_type = config.models["intensity"]
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
        search = pp.load_search(config.output / "msms" / spectra_file.with_suffix(".rescore").name)
        library = pp.merge_spectra_and_peptides(spectra, search)
        annotate_neutral_loss = config.ptm_use_neutral_loss
        aspec = pp.annotate_spectral_library(
            psms=library,
            mass_tol=config.mass_tolerance,
            unit_mass_tol=config.unit_mass_tolerance,
            fragmentation_method=config.fragmentation_method,
            custom_mods=config.unimod_to_mass(),
            annotate_neutral_loss=annotate_neutral_loss,
        )
        aspec.write_as_hdf5(hdf5_path)  # write_metadata_annotation

    return aspec


def _get_best_ce(library: Spectra, spectra_file: Path, config: Config):
    results_dir = config.output / "results"
    results_dir.mkdir(exist_ok=True)
    if config.do_refinement_learning:
        # don't do CE calibration
        return
    if library.obs["FRAGMENTATION"].str.endswith("HCD").any():
        use_ransac_model = config.use_ransac_model
        predictor = pr.Predictor.from_config(config, model_type="intensity")
        alignment_library = predictor.ce_calibration(
            library,
            config.ce_range,
            use_ransac_model,
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

            delta_ce = ransac.predict(library.obs[["MASS", "PRECURSOR_CHARGE"]])
            library.obs["COLLISION_ENERGY"] = np.maximum(0, library.obs["COLLISION_ENERGY"] + delta_ce)

        else:
            ce_alignment = alignment_library.obs.groupby(by=["COLLISION_ENERGY"])["SPECTRAL_ANGLE"].mean()

            best_ce = ce_alignment.idxmax()
            pl.plot_mean_sa_ce(
                sa_ce_df=ce_alignment.to_frame().reset_index(),
                filename=results_dir / f"{spectra_file.stem}_mean_spectral_angle_ce.svg",
            )
            pl.plot_violin_sa_ce(
                sa_ce_df=alignment_library.obs[["COLLISION_ENERGY", "SPECTRAL_ANGLE"]],
                filename=results_dir / f"{spectra_file.stem}_violin_spectral_angle_ce.svg",
            )
            library.obs["COLLISION_ENERGY"] = best_ce
            with open(results_dir / f"{spectra_file.stem}_ce.txt", "w") as f:
                f.write(str(best_ce))
                f.close()
    else:
        best_ce = 35
        library.obs["COLLISION_ENERGY"] = best_ce

        with open(results_dir / f"{spectra_file.stem}_ce.txt", "w") as f:
            f.write(str(best_ce))
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

    elif library_input_type == "internal":
        pass
    else:
        raise ValueError(
            f'Library input type {library_input_type} not understood. Can only be "fasta", "peptides", or "internal".'
        )
    spec_library = pp.gen_lib(library_file)

    pp_and_filter_step = ProcessStep(config.output, "speclib_filtered")

    data_dir = config.output / "data"
    if not pp_and_filter_step.is_done():
        data_dir.mkdir(exist_ok=True)
        spec_library = pp.process_and_filter_spectra_data(
            library=spec_library, model=config.models["intensity"], tmt_label=config.tag
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
                f"Found existing spectral library {out_file}. "
                "Attempting to append missing batches from previous run..."
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
            f"Prediction for {n_failed} / {i+1} batches failed. Check the log to find out why. "
            "Then rerun without changing the config file to append only the missing batches to your output file."
        )
        with open(failed_batch_file, "wb") as fh:
            pickle.dump(failed_batches, fh)
        sys.exit(1)


def generate_spectral_lib(config_path: Union[str, Path]):
    """
    Create a SpectralLibrary object and generate the spectral library.

    # TODO full description
    :param config_path: path to config file
    """
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
            out_file, failed_batch_file, spec_library.obs, batch_size, config.models["intensity"]
        )
        speclib = writer(out_file, mode=mode, min_intensity_threshold=config.min_intensity)
        n_batches = len(batches)

        with Manager() as manager:
            # setup
            shared_queue = manager.Queue(maxsize=config.num_threads)
            prediction_progress = manager.Value("i", 0)
            prediction_failure_progress = manager.Value("i", 0)
            writing_progress = manager.Value("i", 0)

            lock = manager.Lock()
            lock_failure = manager.Lock()

            # Create a pool for producer processes
            predictor_pool = pool.Pool(config.num_threads)

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
    aspec = _annotate_and_get_library(spectra_file, config, tims_meta_file=tims_meta_file)
    if not config.do_refinement_learning:
        _get_best_ce(aspec, spectra_file, config)

    aspec.write_as_hdf5(config.output / "data" / spectra_file.with_suffix(".mzml.pred.hdf5").name)

    ce_calib_step.mark_done()

    return aspec


def run_ce_calibration(
    config_path: Union[str, Path],
):
    """
    Create a CeCalibration object and run the CE calibration.

    # TODO full description
    :param config_path: path to config file
    """
    config = Config()
    config.read(config_path)
    config.check()

    # load spectra file names
    spectra_files = pp.list_spectra(input_dir=config.spectra, input_format=config.spectra_type)

    proc_dir = config.output / "proc"
    proc_dir.mkdir(parents=True, exist_ok=True)

    spectra_files = _preprocess(spectra_files, config)

    if config.num_threads > 1:
        processing_pool = JobPool(processes=config.num_threads)
        for spectra_file in spectra_files:
            processing_pool.apply_async(_ce_calib, [spectra_file, config])
        processing_pool.check_pool()
    else:
        for spectra_file in spectra_files:
            _ce_calib(spectra_file, config)


def _refinement_learn(spectra_files: list[Path], config: Config):
    refinement_step = ProcessStep(config.output, "refinement_learning")
    if refinement_step.is_done():
        return

    libraries = [_ce_calib(spectra_file, config) for spectra_file in spectra_files]

    if config.download_baseline_intensity_predictor:
        baseline_model_path = config.output / "data/dlomix/prosit_baseline_model.keras"
        download_new_baseline_model = True
    else:
        baseline_model_path = Path(config.models["intensity"])
        download_new_baseline_model = False

    pr.dlomix.refine_intensity_predictor(
        baseline_model_path=baseline_model_path,
        libraries=libraries,
        config=config,
        data_directory=config.output / "data/dlomix",
        result_directory=config.output / "results/dlomix",
        dataset_name="refinement_dataset",
        model_name="refined",
        download_new_baseline_model=download_new_baseline_model,
    )

    refinement_step.mark_done()


def _calculate_features(spectra_file: Path, config: Config):
    library = _ce_calib(spectra_file, config)
    calc_feature_step = ProcessStep(config.output, "calculate_features." + spectra_file.stem)
    if calc_feature_step.is_done():
        return

    predict_step = ProcessStep(config.output, "predict." + spectra_file.stem)
    if not predict_step.is_done():

        if "alphapept" in config.models["intensity"].lower():
            chunk_idx = list(group_iterator(df=library.obs, group_by_column="PEPTIDE_LENGTH"))
        else:
            chunk_idx = None
        if config.do_refinement_learning:
            intensity_predictor = pr.Predictor.from_dlomix(
                model_type="intensity",
                model_path=config.output / "data/dlomix/refined.keras",
                output_path=config.output / "data/dlomix/",
                batch_size=config.dlomix_inference_batch_size,
            )
        else:
            intensity_predictor = pr.Predictor.from_config(config, model_type="intensity")

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
        regression_method=config.curve_fitting_method,
    )
    re.generate_features(
        library=library,
        search_type="rescore",
        output_file=fdr_dir / spectra_file.with_suffix(".rescore.tab").name,
        additional_columns=config.use_feature_cols,
        all_features=config.all_features,
        regression_method=config.curve_fitting_method,
        add_neutral_loss_features=add_neutral_loss_features,
        remove_miss_cleavage_features=remove_miss_cleavage_features,
    )

    calc_feature_step.mark_done()


def _rescore(fdr_dir: Path, config: Config):
    """
    High level rescore function for original and rescore.

    :param fdr_dir: the output directory
    :param config: the configuration object
    :raises ValueError: if the provided fdr estimation method in the config is not recognized
    """
    rescore_original_step = ProcessStep(config.output, f"{config.fdr_estimation_method}_original")
    rescore_prosit_step = ProcessStep(config.output, f"{config.fdr_estimation_method}_prosit")

    if config.fdr_estimation_method == "percolator":
        if not rescore_original_step.is_done():
            re.rescore_with_percolator(input_file=fdr_dir / "original.tab", output_folder=fdr_dir)
            rescore_original_step.mark_done()
        if not rescore_prosit_step.is_done():
            logger.info("Start percolator rescoring")
            logger.info(config.ptm_localization)
            if config.ptm_localization:
                _ptm_localization_rescore(fdr_dir, config)
            else:
                re.rescore_with_percolator(input_file=fdr_dir / "rescore.tab", output_folder=fdr_dir)
                rescore_prosit_step.mark_done()
    elif config.fdr_estimation_method == "mokapot":
        if not rescore_original_step.is_done():
            re.rescore_with_mokapot(input_file=fdr_dir / "original.tab", output_folder=fdr_dir)
            rescore_original_step.mark_done()
        if not rescore_prosit_step.is_done():
            re.rescore_with_mokapot(input_file=fdr_dir / "rescore.tab", output_folder=fdr_dir)
            rescore_prosit_step.mark_done()
    else:
        raise ValueError(
            'f{config.fdr_estimation_method} is not a valid rescoring tool, use either "percolator" or "mokapot"'
        )


def _ptm_localization_rescore(fdr_dir: Path, config: Config):
    """
     Helper function for running percolator to do PTM localization.

    :param fdr_dir: the output directory
    :param config: the configuration object
    """
    df_rescore = pd.read_csv(fdr_dir / "rescore.tab", sep="\t")
    unimod_id = config.ptm_unimod_id
    df_rescore["id"] = df_rescore["filename"] + df_rescore["ScanNr"].astype(str)
    df_rescore_fil = df_rescore[df_rescore["Peptide"].str.contains("UNIMOD:" + str(unimod_id))]
    df_rescore = df_rescore[df_rescore["id"].isin(df_rescore_fil["id"].tolist())]
    df_rescore.drop(columns=["id"], inplace=True)
    df_rescore.to_csv(fdr_dir / "rescore.tab", sep="\t", index=False)
    new_rescore_dir = fdr_dir / "localize_mod"
    new_rescore_dir.mkdir(parents=True, exist_ok=True)

    if unimod_id == 7:
        re.rescore_with_percolator(input_file=fdr_dir / "rescore.tab", output_folder=fdr_dir)
        df_rescore_psms_targets = pd.read_csv(fdr_dir / "rescore.percolator.psms.txt", sep="\t")
        df_rescore_psms_decoys = pd.read_csv(fdr_dir / "rescore.percolator.decoy.psms.txt", sep="\t")
        df_rescore_psms_targets = df_rescore_psms_targets[
            df_rescore_psms_targets["peptide"].apply(lambda x: "UNIMOD:" + str(unimod_id) in x)
        ]
        df_rescore_psms_decoys = df_rescore_psms_decoys[
            df_rescore_psms_decoys["peptide"].apply(lambda x: "UNIMOD:" + str(unimod_id) in x)
        ]
        df_rescore_psms = pd.concat([df_rescore_psms_targets, df_rescore_psms_decoys])
        df_rescore_psms = df_rescore_psms[["PSMId"]]
        df_rescore_psms.rename(columns={"PSMId": "SpecId"}, inplace=True)
        df_rescore = df_rescore.merge(df_rescore_psms, on="SpecId", how="inner")
        df_rescore.to_csv(new_rescore_dir / "rescore.tab", sep="\t", index=False)
        re.rescore_with_percolator(input_file=new_rescore_dir / "rescore.tab", output_folder=new_rescore_dir)
    else:
        re.rescore_with_percolator(input_file=fdr_dir / "rescore.tab", output_folder=new_rescore_dir)


def run_rescoring(config_path: Union[str, Path]):
    """
    Create a ReScore object and run the rescoring.

    # TODO full description
    :param config_path: path to config file
    """
    config = Config()
    config.read(config_path)
    config.check()

    # load spectra file names
    spectra_files = pp.list_spectra(input_dir=config.spectra, input_format=config.spectra_type)

    proc_dir = config.output / "proc"
    proc_dir.mkdir(parents=True, exist_ok=True)

    spectra_files = _preprocess(spectra_files, config)

    # TODO is this the most elegant way to multi-thread CE calibration before running refinement learning?
    # Should we store the returned libraries and pass them to _calculate_features and _refinement_learn instead of
    # _ce_calib returning cached outputs?
    if config.num_threads > 1:
        processing_pool = JobPool(processes=config.num_threads)
        for spectra_file in spectra_files:
            _ = processing_pool.apply_async(_ce_calib, [spectra_file, config])
        processing_pool.check_pool()
    else:
        for spectra_file in spectra_files:
            _ = _ce_calib(spectra_file, config)

    if config.do_refinement_learning:
        _refinement_learn(spectra_files, config)

    if config.num_threads > 1:
        processing_pool = JobPool(processes=config.num_threads)
        for spectra_file in spectra_files:
            processing_pool.apply_async(_calculate_features, [spectra_file, config])
        processing_pool.check_pool()
    else:
        for spectra_file in spectra_files:
            _calculate_features(spectra_file, config)

    # prepare rescoring

    fdr_dir = config.output / "results" / config.fdr_estimation_method

    original_tab_files = [fdr_dir / spectra_file.with_suffix(".original.tab").name for spectra_file in spectra_files]
    rescore_tab_files = [fdr_dir / spectra_file.with_suffix(".rescore.tab").name for spectra_file in spectra_files]

    prepare_tab_original_step = ProcessStep(config.output, f"{config.fdr_estimation_method}_prepare_tab_original")
    prepare_tab_rescore_step = ProcessStep(config.output, f"{config.fdr_estimation_method}_prepare_tab_prosit")

    if not prepare_tab_original_step.is_done():
        logger.info("Merging input tab files for rescoring without peptide property prediction")
        re.merge_input(tab_files=original_tab_files, output_file=fdr_dir / "original.tab")
        prepare_tab_original_step.mark_done()

    if not prepare_tab_rescore_step.is_done():
        logger.info("Merging input tab files for rescoring with peptide property prediction")
        re.merge_input(tab_files=rescore_tab_files, output_file=fdr_dir / "rescore.tab")
        prepare_tab_rescore_step.mark_done()

    # rescoring
    _rescore(fdr_dir, config)

    # plotting
    logger.info("Generating summary plots...")
    if not config.ptm_localization:
        pl.plot_all(fdr_dir)
    logger.info("Finished rescoring.")

    if config.quantification:
        logger.info("Starting quantification")
        # method contains picked-group-FDR call
        apply_quant(config)
        logger.info("Finished quantification")


def run_job(config_path: Union[str, Path]):
    """
    Run oktoberfest based on job type given in the config file.

    :param config_path: Path to config file
    :raises ValueError: In case the job_type in the provided config file is not known
    """
    conf = Config()
    conf.read(config_path)
    conf.check()

    output_folder = conf.output
    job_type = conf.job_type

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
    logger.info(json.dumps(conf.data, indent=4))

    try:
        if job_type == "SpectralLibraryGeneration":
            generate_spectral_lib(config_path)
        elif job_type == "CollisionEnergyCalibration":
            run_ce_calibration(config_path)
        elif job_type == "Rescoring":
            run_rescoring(config_path)
        else:
            raise ValueError(f"Unknown job_type in config: {job_type}")
    finally:
        file_handler.close()
        base_logger.removeHandler(file_handler)
