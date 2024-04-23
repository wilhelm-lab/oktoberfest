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
from typing import Dict, List, Optional, Tuple, Type, Union
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression, RANSACRegressor
from spectrum_io.spectral_library import MSP, DLib, SpectralLibrary, Spectronaut
from tqdm.auto import tqdm

from oktoberfest import __copyright__, __version__
from oktoberfest import plotting as pl
from oktoberfest import predict as pr
from oktoberfest import preprocessing as pp
from oktoberfest import rescore as re

from .data.spectra import FragmentType, Spectra
from .utils import Config, JobPool, ProcessStep

logger = logging.getLogger(__name__)


def _make_predictions_error_callback(failure_progress_tracker, failure_lock, error):
    logger.error(
        f"Prediction failed due to: {error} Batch will be missing from output file. "
        "DO NOT STOP THIS RUN: The index of the batch is stored and your output file will be appended "
        "by the missing batches if you rerun without changing your config file after this run is completed."
    )
    with failure_lock:
        failure_progress_tracker.value += 1


def _make_predictions(int_model, irt_model, predict_kwargs, queue_out, progress, lock, batch_df):
    predictions = {
        **pr.predict(batch_df, model_name=int_model, **predict_kwargs),
        **pr.predict(batch_df, model_name=irt_model, **predict_kwargs),
    }
    queue_out.put((predictions, batch_df))
    with lock:
        progress.value += 1


def _preprocess(spectra_files: List[Path], config: Config) -> List[Path]:
    preprocess_search_step = ProcessStep(config.output, "preprocessing_search")
    if not preprocess_search_step.is_done():
        # load search results
        if not config.search_results_type == "internal":
            logger.info(f"Converting search results from {config.search_results} to internal search result.")

            msms_output = config.output / "msms"
            msms_output.mkdir(exist_ok=True)
            internal_search_file = msms_output / "msms.prosit"
            tmt_label = config.tag

            search_results = pp.convert_search(
                input_path=config.search_results,
                search_engine=config.search_results_type,
                tmt_label=tmt_label,
                output_file=internal_search_file,
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
        search_results = pp.filter_peptides_for_model(peptides=search_results, model=config.models["intensity"])

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
        search = pp.load_search(config.output / "msms" / spectra_file.with_suffix(".rescore").name)
        library = pp.merge_spectra_and_peptides(spectra, search)
        if "xl" in config.models["intensity"].lower():
            aspec = pp.annotate_spectral_library_xl(library, mass_tol=config.mass_tolerance, unit_mass_tol=config.unit_mass_tolerance)
        else:
            paspecp.annotate_spectral_library(library, mass_tol=config.mass_tolerance, unit_mass_tol=config.unit_mass_tolerance)
        aspec.write_as_hdf5(hdf5_path)  # write_metadata_annotation

    return aspec


def _get_best_ce(library: Spectra, spectra_file: Path, config: Config):
    results_dir = config.output / "results"
    results_dir.mkdir(exist_ok=True)
    if (library.obs["FRAGMENTATION"] == "HCD").any():
        server_kwargs = {
            "server_url": config.prediction_server,
            "ssl": config.ssl,
            "model_name": config.models["intensity"],
        }
        use_ransac_model = config.use_ransac_model
        if "xl" in config.models["intensity"].lower():
            alignment_library = pr.ce_calibration(library, config.ce_range, use_ransac_model, xl =True, **server_kwargs )
        else:
            alignment_library = pr.ce_calibration(library, config.ce_range, use_ransac_model, **server_kwargs)

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
    if library_input_type == "fasta":
        digest_step = ProcessStep(config.output, "speclib_digested")
        library_file = config.output / "prosit_input.csv"
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
            metadata = pp.generate_metadata(
                peptides=list(peptide_dict.keys()),
                collision_energy=config.collision_energy,
                precursor_charge=config.precursor_charge,
                fragmentation=config.fragmentation,
                proteins=list(peptide_dict.values()),
            )
            library_file = config.output / "prosit_input.csv"
            metadata.to_csv(library_file, sep=",", index=None)
            digest_step.mark_done()
    elif library_input_type == "peptides":
        library_file = config.library_input
    else:
        raise ValueError(f'Library input type {library_input_type} not understood. Can only be "fasta" or "peptides".')
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


def _get_writer_and_output(results_path: Path, output_format: str) -> Tuple[Type[SpectralLibrary], Path]:
    if output_format == "msp":
        return MSP, results_path / "myPrositLib.msp"
    elif output_format == "spectronaut":
        return Spectronaut, results_path / "myPrositLib.csv"
    elif output_format == "dlib":
        return DLib, results_path / "myPrositLib.dlib"
    else:
        raise ValueError(f"{output_format} is not supported as spectral library type")


def _get_batches_and_mode(out_file: Path, failed_batch_file: Path, no_of_spectra: int, batchsize: int):
    if out_file.is_file():
        if failed_batch_file.is_file():
            with open(failed_batch_file, "rb") as fh:
                batches = pickle.load(fh)
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
        batches = range(ceil(no_of_spectra / batchsize))
        mode = "w"

    return batches, mode


def _update(pbar: tqdm, postfix_values: Dict[str, int]):
    total_val = sum(postfix_values.values())
    if total_val > pbar.n:
        pbar.set_postfix(**postfix_values)
        pbar.n = total_val
        pbar.refresh()


def _check_write_failed_batch_file(failed_batch_file: Path, n_failed: int, results: List[pool.AsyncResult]):
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

    spec_library = _speclib_from_digestion(config)

    server_kwargs = {
        "server_url": config.prediction_server,
        "ssl": config.ssl,
        "disable_progress_bar": True,
    }

    speclib_written_step = ProcessStep(config.output, "speclib_written")
    if not speclib_written_step.is_done():
        results_path = config.output / "results"
        results_path.mkdir(exist_ok=True)

        batchsize = config.batch_size
        failed_batch_file = config.output / "data" / "speclib_failed_batches.pkl"
        writer, out_file = _get_writer_and_output(results_path, config.output_format)
        batches, mode = _get_batches_and_mode(out_file, failed_batch_file, spec_library.n_obs, batchsize)
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
                args=(
                    shared_queue,
                    writing_progress,
                ),
            )

            try:
                results = []
                for i in batches:
                    result = predictor_pool.apply_async(
                        _make_predictions,
                        (
                            config.models["intensity"],
                            config.models["irt"],
                            server_kwargs,
                            shared_queue,
                            prediction_progress,
                            lock,
                            spec_library.obs[i * batchsize : (i + 1) * batchsize],
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


def _calculate_features(spectra_file: Path, config: Config, xl: bool = False):
    library = _ce_calib(spectra_file, config)

    calc_feature_step = ProcessStep(config.output, "calculate_features." + spectra_file.stem)
    if calc_feature_step.is_done():
        return

    predict_step = ProcessStep(config.output, "predict." + spectra_file.stem)
    if not predict_step.is_done():

        predict_kwargs = {
            "server_url": config.prediction_server,
            "ssl": config.ssl,
        }
        if xl:
            pred_intensities_a, pred_intensities_b  = pr.predict(
            data=library.obs,
            model_name=config.models["intensity"],
            xl=True,
            **predict_kwargs,
        )
            library.add_matrix(pred_intensities_a["intensities"], FragmentType.PRED_A)
            library.add_matrix(pred_intensities_b["intensities"], FragmentType.PRED_B)
            library.write_as_hdf5(config.output / "data" / spectra_file.with_suffix(".mzml.pred.hdf5").name)
            predict_step.mark_done()
        else:
            pred_intensities = pr.predict(
            data=library.obs,
            model_name=config.models["intensity"],
            **predict_kwargs,
        )
            pred_irts = pr.predict(data=library.obs, model_name=config.models["irt"], **predict_kwargs)
            library.add_matrix(pred_intensities["intensities"], FragmentType.PRED)
            library.add_column(pred_irts["irt"].squeeze(), name="PREDICTED_IRT")
            library.write_as_hdf5(config.output / "data" / spectra_file.with_suffix(".mzml.pred.hdf5").name)
            predict_step.mark_done()

    # produce percolator tab files
    fdr_dir = config.output / "results" / config.fdr_estimation_method
    fdr_dir.mkdir(exist_ok=True)

    re.generate_features(
        library=library,
        search_type="original",
        output_file=fdr_dir / spectra_file.with_suffix(".original.tab").name,
        all_features=config.all_features,
        xl=xl,
        regression_method=config.curve_fitting_method,
    )
    re.generate_features(
        library=library,
        search_type="rescore",
        output_file=fdr_dir / spectra_file.with_suffix(".rescore.tab").name,
        all_features=config.all_features,
        xl=xl,
        regression_method=config.curve_fitting_method,
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

def xl_fdr(df: pd.DataFrame, score: str) -> pd.DataFrame: 
    """
    "calculate and add fdr_xl to the DataFrame : (TD-DD)/(TT)".  

    :param df: DataFrame containing the data.
    :param score: Column name containing the scores used for calculating FDR 
    :return: DataFrame with the new column 'fdr'
    """
    df = df.sort_values(by=score, ascending=False)
    df['TD_sum'] = (df['label'] == 'TD').cumsum()
    df['DD_sum'] = (df['label'] == 'DD').cumsum()
    df['TT_sum'] = (df['label'] == 'TT').cumsum()
    df['mokapot q-value'] = (df['TD_sum'] - df['DD_sum']) / df['TT_sum']
    df.loc[df['TT_sum'] == 0, 'mokapot q-value'] = 0  # Handling division by zero
    df = df.drop(['TD_sum', 'DD_sum', 'TT_sum'], axis=1)
    return df


def xl_preprocessing_plot_unique_xl(featrures_dir: Path, mokapot_csms: pd.DataFrame, original_or_rescore: str):
    columns_to_keep_mokapot = ["SpecId_a", "label", "ScanNr_a", "filename_a", "Peptide_a", "mokapot score", "mokapot q-value",  "mokapot q-value_a", "mokapot q-value_b", "Proteins_b", "mod_pep_a_b", "mod_pep_b_b"]
    mokapot_csms = mokapot_csms[columns_to_keep_mokapot]
    mokapot_csms = mokapot_csms.rename(columns={'SpecId_a': 'SpecId', 'ScanNr_a': 'ScanNr', 'filename_a': 'filename', 'Peptide_a': 'Peptide', 'Proteins_b': 'Proteins', 'mod_pep_a_b': 'mod_pep_a', 'mod_pep_b_b': 'mod_pep_b'})
    mokapot_csms['label'] = mokapot_csms['label'].replace({'TT': True, 'TD': False, 'DD': False})
    mokapot_csms['base_pep_a'] = mokapot_csms['mod_pep_a'].str.replace(r'\[.*?\]', '', regex=True)
    mokapot_csms['base_pep_b'] = mokapot_csms['mod_pep_b'].str.replace(r'\[.*?\]', '', regex=True)
    mokapot_csms.to_csv("/cmnfs/data/proteomics/XL/Ribosome/mokapot_csms_10.csv")
    mokapot_csms = mokapot_csms.sort_values(by='mokapot score', ascending = False)
    swap_mask = mokapot_csms['base_pep_a'] > mokapot_csms['base_pep_b']
    mokapot_csms.loc[swap_mask, 'min_seq'] = mokapot_csms['base_pep_b']
    mokapot_csms.loc[swap_mask, 'max_seq'] = mokapot_csms['base_pep_a']
    mokapot_csms.loc[~swap_mask, 'min_seq'] = mokapot_csms['base_pep_a']
    mokapot_csms.loc[~swap_mask, 'max_seq'] = mokapot_csms['base_pep_b']
    mokapot_uniq_xls = mokapot_csms.drop_duplicates(subset=['min_seq', 'max_seq'])
    xl_fdr(mokapot_uniq_xls, "mokapot score")
    mokapot_uniq_xls_true = mokapot_uniq_xls[mokapot_uniq_xls['label'] == True]
    mokapot_uniq_xls_decoy = mokapot_uniq_xls[mokapot_uniq_xls['label'] == False]
    if original_or_rescore == "original":
        mokapot_uniq_xls_true.to_csv(featrures_dir + "/original.mokapot.peptides.txt", sep='\t', index=False)
        mokapot_uniq_xls_decoy.to_csv(featrures_dir + "/original.mokapot.decoy.peptides.txt", sep='\t', index=False)
    else:
        mokapot_uniq_xls_true.to_csv(featrures_dir + "/rescore.mokapot.peptides.txt", sep='\t', index=False)
        mokapot_uniq_xls_decoy.to_csv(featrures_dir + "/rescore.mokapot.decoy.peptides.txt", sep='\t', index=False)

    
def xl_preprocessing_plot_csm(featrures_dir: Path, mokapot_csms: pd.DataFrame, original_or_rescore: str):
    columns_to_keep_mokapot = ["SpecId_a", "label", "ScanNr_a", "filename_a", "Peptide_a", "mokapot score", "mokapot q-value",  "mokapot q-value_a", "mokapot q-value_b", "Proteins_b"]
    mokapot_csms = mokapot_csms[columns_to_keep_mokapot]
    mokapot_csms = mokapot_csms.rename(columns={'SpecId_a': 'SpecId', 'ScanNr_a': 'ScanNr', 'filename_a': 'filename', 'Peptide_a': 'Peptide', 'Proteins_b': 'Proteins'})
    mokapot_csms['label'] = mokapot_csms['label'].replace({'TT': True, 'TD': False, 'DD': False})
    mokapot_csms_true = mokapot_csms[mokapot_csms['label'] == True]
    mokapot_csms_decoy = mokapot_csms[mokapot_csms['label'] == False]
    if original_or_rescore == "original":
        mokapot_csms_true.to_csv(featrures_dir + "/original.mokapot.psms.txt", sep='\t', index=False)
        mokapot_csms_decoy.to_csv(featrures_dir + "/original.mokapot.decoy.psms.txt", sep='\t', index=False)
    else:
        mokapot_csms_true.to_csv(featrures_dir + "/rescore.mokapot.psms.txt", sep='\t', index=False)
        mokapot_csms_decoy.to_csv(featrures_dir + "/rescore.mokapot.decoy.psms.txt", sep='\t', index=False)
    

def xl_psm_to_csm(featrures_dir: Path, original_or_rescore: str ):  
    def get_label(row):
        if row['is_decoy_p1_a'] == 'False' and row['is_decoy_p2_a'] == 'False':
            return 'TT'
        elif row['is_decoy_p1_a'] == 'True' and row['is_decoy_p2_a'] == 'True':
            return 'DD'
        else:
            return 'TD'
    
    def min_score(row):
        return min(row['mokapot score_a'], row['mokapot score_b'])
    
    if original_or_rescore == "original":
        decoy_psms = pd.read_csv(featrures_dir + "/original.mokapot.decoy.psms.txt", delimiter='\t') 
        target_psms = pd.read_csv(featrures_dir + "/original.mokapot.psms.txt", delimiter='\t') 
    else:
        decoy_psms = pd.read_csv(featrures_dir + "/rescore.mokapot.decoy.psms.txt", delimiter='\t') 
        target_psms = pd.read_csv(featrures_dir + "/rescore.mokapot.psms.txt", delimiter='\t') 

    split_data = target_psms['SpecId'].str.rsplit('-', n=10, expand=True)
    new_columns = ["raw_file","scan_number", "mod_pep_a", "mod_pep_b", "charge",  
               "decoy_p1", "is_decoy_p1", "decoy_p2", "is_decoy_p2", "index"]
    split_data.columns = new_columns
    df_psm_target = pd.concat([target_psms, split_data], axis=1)
    split_data = decoy_psms['SpecId'].str.rsplit('-', n=10, expand=True)
    new_columns = ["raw_file","scan_number", "mod_pep_a", "mod_pep_b", "charge", 
                "decoy_p1", "is_decoy_p1", "decoy_p2", "is_decoy_p2", "index"]
    split_data.columns = new_columns
    df_psm_decoy = pd.concat([decoy_psms, split_data], axis=1)
    df_mokapot = pd.concat([df_psm_decoy , df_psm_target], axis = 0)
    df_mokapot[['index_csm', '_', 'which_pep',]] = df_mokapot['index'].str.split('_', expand=True)
    df_mokapot.drop(columns=['index', '_', 'decoy_p1', 'decoy_p2'], inplace=True)
    df_pep_1 = df_mokapot[df_mokapot['which_pep'] == "1"].copy()
    df_pep_2 = df_mokapot[df_mokapot['which_pep'] == "2"].copy()
    df_pep_1.drop(columns=['which_pep'], inplace=True)
    df_pep_2.drop(columns=['which_pep'], inplace=True)
    df_pep_1.columns = [col + '_a' if col != 'index_csm' else col for col in df_pep_1.columns]
    df_pep_2.columns = [col + '_b' if col != 'index_csm' else col for col in df_pep_2.columns]
    mokapot_csms = pd.merge(df_pep_1, df_pep_2, on='index_csm')
    mokapot_csms['mokapot score'] = mokapot_csms.apply(min_score, axis=1)
    mokapot_csms['label'] = mokapot_csms.apply(get_label, axis=1)
    mokapot_csms.to_csv("/cmnfs/data/proteomics/XL/Ribosome/mokapot_csms.csv")
    return mokapot_csms


def prepare_rescore_xl_psm_level(featrures_dir: Path, original_or_rescore: str ):
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
   
    columns_to_remove_psm_rescore= ["run_name",
            "scan_number",
            "precursor_mass",
            "precursor_charge",
            "crosslinker_name",
            "base_sequence_p1",
            "aa_len_p1",
            "link_pos_p1",
            "linked_aa_p1",
            "mods_p1",
            "mod_pos_p1",
            "base_sequence_p2",
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
            "Unnamed: 0"
    ]
    columns_to_remove_psm_original = list(set(columns_to_remove_psm_rescore + ["match_score"]) - set(["spectral_angle"]))

    if original_or_rescore == "original":
        columns_to_remove_psm = columns_to_remove_psm_original
        rescore_tab_file =  pd.read_csv(featrures_dir + "/original.tab", sep="\t")
    else:
        columns_to_remove_psm = columns_to_remove_psm_rescore
        rescore_tab_file =  pd.read_csv(featrures_dir + "/rescore.tab", sep="\t")


    rescore_tab_file.drop(columns=columns_to_remove_psm, inplace=True)
    rescore_tab_file = rescore_tab_file.fillna(0)
    string_columns = ['SpecId', 'Peptide', 'Proteins', 'protein_p1', 'protein_p2', 'filename']  
    rescore_tab_file_numeric_columns = [col for col in  rescore_tab_file.columns if col not in string_columns]
    rescore_tab_file[rescore_tab_file_numeric_columns] = rescore_tab_file[rescore_tab_file_numeric_columns].apply(pd.to_numeric, errors='coerce')
    rescore_tab_file = rescore_tab_file.reset_index(drop=True)
    rescore_tab_file['Proteins'] = 'p1_' + rescore_tab_file['protein_p1'].astype(str) + '_p2_' + rescore_tab_file['protein_p2'].astype(str)
    rescore_tab_file['SpecId'] = (
    rescore_tab_file['SpecId'] + '-decoy_p1-' +
    rescore_tab_file['decoy_p1'].astype(str) + '-decoy_p2-' +
    rescore_tab_file['decoy_p2'].astype(str) + '-' +
    rescore_tab_file.index.astype(str)         # Adding index to SpecId
)
    rescore_tab_file.drop(columns=[ 'protein_p1', 'protein_p2', 'decoy_p1', 'decoy_p2'], inplace=True)
    columns_feature_pep_a = [col for col in rescore_tab_file.columns if col.endswith('_a') or col.endswith('_A') or not (col.endswith('_b') or col.endswith('_B'))]
    columns_feature_pep_b = [col for col in rescore_tab_file.columns if col.endswith('_b') or col.endswith('_B') or not (col.endswith('_a') or col.endswith('_A'))]
    rescore_tab_file_a = rescore_tab_file.loc[:, columns_feature_pep_a]
    rescore_tab_file_b = rescore_tab_file.loc[:, columns_feature_pep_b]
    rescore_tab_file_a['SpecId'] = rescore_tab_file_a['SpecId'] + '_' + rescore_tab_file_a.index.astype(str)
    rescore_tab_file_b['SpecId'] = rescore_tab_file_b['SpecId'] + '_' + rescore_tab_file_b.index.astype(str)
    #1 means pep a and 2 means pep b
    rescore_tab_file_a['SpecId'] = rescore_tab_file_a['SpecId'] + '_1' 
    rescore_tab_file_b['SpecId'] = rescore_tab_file_b['SpecId'] + '_2' 
    rescore_tab_file_a.columns = [col[:-2] if col.endswith(('_a', '_A')) else col for col in rescore_tab_file_a.columns]
    rescore_tab_file_b.columns = [col[:-2] if col.endswith(('_b', '_B')) else col for col in rescore_tab_file_b.columns]
    rescore_tab_file_a['Label'] = rescore_tab_file_a['SpecId'].apply(extract_label_pep_a)
    rescore_tab_file_b['Label'] = rescore_tab_file_b['SpecId'].apply(extract_label_pep_b)
    string_columns = ['SpecId', 'Peptide', 'Proteins', 'filename']  
    rescore_tab_file_numeric_columns = [col for col in  rescore_tab_file_a.columns if col not in string_columns]
    rescore_tab_file_a[rescore_tab_file_numeric_columns] = rescore_tab_file_a[rescore_tab_file_numeric_columns].apply(pd.to_numeric, errors='coerce')
    rescore_tab_file_b[rescore_tab_file_numeric_columns] = rescore_tab_file_b[rescore_tab_file_numeric_columns].apply(pd.to_numeric, errors='coerce')
    # change ExpMass of rescore_tab_file_a
    max_ExpMass = rescore_tab_file_a['ExpMass'].max()
    rescore_tab_file_b['ExpMass'] += max_ExpMass
    input_mokapot_psm = pd.concat([rescore_tab_file_a, rescore_tab_file_b], axis=0, ignore_index=True)
    input_mokapot_psm['Proteins'].fillna("unknown", inplace=True)
    return input_mokapot_psm


def run_rescoring(config_path: Union[str, Path]):
    """
    Create a ReScore object and run the rescoring.

    # TODO full description
    :param config_path: path to config file
    """
    config = Config()
    config.read(config_path)

    # load spectra file names
    spectra_files = pp.list_spectra(input_dir=config.spectra, input_format=config.spectra_type)

    proc_dir = config.output / "proc"
    proc_dir.mkdir(parents=True, exist_ok=True)

    spectra_files = _preprocess(spectra_files, config)

    if config.num_threads > 1:
        processing_pool = JobPool(processes=config.num_threads)
        for spectra_file in spectra_files:
            if "xl" in config.models["intensity"].lower():
                processing_pool.apply_async(_calculate_features, [spectra_file, config], xl=True)
            else:
                processing_pool.apply_async(_calculate_features, [spectra_file, config])
        processing_pool.check_pool()
    else:
        for spectra_file in spectra_files:
            if "xl" in config.models["intensity"].lower():
                _calculate_features(spectra_file, config, xl=True)
            else:
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
    if "xl" in config.models["intensity"].lower(): #xl-psm-level
        input_mokapot_psm_rescore = prepare_rescore_xl_psm_level(str(fdr_dir), "rescore")
        input_mokapot_psm_rescore.to_csv(str(fdr_dir) + "/rescore.tab", sep="\t", index=None)
        input_mokapot_psm_original = prepare_rescore_xl_psm_level(str(fdr_dir), "original")
        input_mokapot_psm_original.to_csv(str(fdr_dir) + "/original.tab", sep="\t", index=None)
        _rescore(fdr_dir, config)
        mokapot_csms_rescore = xl_psm_to_csm(str(fdr_dir), "rescore")
        mokapot_csms_original = xl_psm_to_csm(str(fdr_dir), "original")
        mokapot_csms_rescore = xl_fdr(mokapot_csms_rescore, score = "mokapot score")
        mokapot_csms_original = xl_fdr(mokapot_csms_original, score = "mokapot score")
        xl_preprocessing_plot_csm(str(fdr_dir), mokapot_csms_rescore, "rescore")
        xl_preprocessing_plot_csm(str(fdr_dir), mokapot_csms_original, "original")
        xl_preprocessing_plot_unique_xl(str(fdr_dir), mokapot_csms_rescore, "rescore")
        xl_preprocessing_plot_unique_xl(str(fdr_dir), mokapot_csms_original, "original")

    else:
        _rescore(fdr_dir, config)

    # plotting
    logger.info("Generating summary plots...")
    pl.plot_all(fdr_dir)

    logger.info("Finished rescoring.")


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
