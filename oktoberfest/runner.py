import datetime
import json
import logging
import sys
import time
import traceback
from multiprocessing import Manager, Pool, Process
from pathlib import Path
from typing import List, Type, Union

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

from .data.spectra import FragmentType, Spectra
from .utils import Config, JobPool, ProcessStep

logger = logging.getLogger(__name__)


def _make_predictions(batch_df, int_model, irt_model, server_kwargs, queue_out, progress, lock):
    predictions = pr.predict(batch_df, model_name=int_model, disable_progress_bar=True, **server_kwargs)
    predictions |= pr.predict(batch_df, model_name=irt_model, disable_progress_bar=True, **server_kwargs)
    with lock:
        progress.value += 1
    queue_out.put((predictions, batch_df))


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

            pp.convert_search(
                input_path=config.search_results,
                output_file=internal_search_file,
                search_engine=config.search_results_type,
                tmt_label=tmt_label,
            )
        else:
            internal_search_file = config.search_results
        search_results = pp.load_search(internal_search_file)
        logger.info(f"Read {len(search_results)} PSMs from {internal_search_file}")

        # filter search results
        search_results = pp.filter_peptides_for_model(peptides=search_results, model=config.models["intensity"])

        # split search results
        filenames_found = pp.split_search(
            search_results=search_results,
            output_dir=config.output / "msms",
            filenames=[spectra_file.stem for spectra_file in spectra_files],
        )
        preprocess_search_step.mark_done()
    else:
        filenames_found = [msms_file.stem for msms_file in (config.output / "msms").glob("*rescore")]

    spectra_files_to_return = []
    for spectra_file in spectra_files:
        if spectra_file.stem in filenames_found:
            spectra_files_to_return.append(spectra_file)

    return spectra_files_to_return


def _annotate_and_get_library(spectra_file: Path, config: Config) -> Spectra:
    data_dir = config.output / "data"
    data_dir.mkdir(exist_ok=True)
    hdf5_path = data_dir / spectra_file.with_suffix(".mzml.hdf5").name
    if hdf5_path.is_file():
        library = Spectra.from_hdf5(hdf5_path)
    else:
        spectra_dir = config.output / "spectra"
        spectra_dir.mkdir(exist_ok=True)
        format_ = spectra_file.suffix.lower()
        if format_ == ".raw":
            file_to_load = spectra_dir / spectra_file.with_suffix(".mzML").name
            pp.convert_spectra_to_mzml(spectra_file, file_to_load, thermo_exe=config.thermo_exe)
        elif format_ in [".mzml", ".pkl"]:
            file_to_load = spectra_file
        elif format_ == ".d":
            raise NotImplementedError("Bruker file conversion will be available in spectrum-io v0.3.5")
        spectra = pp.load_spectra(file_to_load)
        search = pp.load_search(config.output / "msms" / spectra_file.with_suffix(".rescore").name)
        library = pp.merge_spectra_and_peptides(spectra, search)
        pp.annotate_spectral_library(library, mass_tol=config.mass_tolerance, unit_mass_tol=config.unit_mass_tolerance)
        library.write_as_hdf5(hdf5_path)  # write_metadata_annotation

    return library


def _get_best_ce(library: Spectra, spectra_file: Path, config: Config):
    results_dir = config.output / "results"
    results_dir.mkdir(exist_ok=True)
    if (library.spectra_data["FRAGMENTATION"] == "HCD").any():
        server_kwargs = {
            "server_url": config.prediction_server,
            "ssl": config.ssl,
            "model_name": config.models["intensity"],
        }
        use_ransac_model = config.use_ransac_model
        library.spectra_data.rename(
            columns={
                "MODIFIED_SEQUENCE": "peptide_sequences",
                "PRECURSOR_CHARGE": "precursor_charges",
                "COLLISION_ENERGY": "collision_energies",
                "FRAGMENTATION": "fragmentation_types",
            },
            inplace=True,
        )
        alignment_library = pr.ce_calibration(library, config.ce_range, use_ransac_model, **server_kwargs)
        library.spectra_data.rename(
            columns={
                "peptide_sequences": "MODIFIED_SEQUENCE",
                "precursor_charges": "PRECURSOR_CHARGE",
                "collision_energies": "COLLISION_ENERGY",
                "fragmentation_types": "FRAGMENTATION",
            },
            inplace=True,
        )

        if use_ransac_model:
            logger.info("Performing RANSAC regression")
            calib_group = (
                alignment_library.spectra_data.groupby(
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

            delta_ce = ransac.predict(library.spectra_data[["MASS", "PRECURSOR_CHARGE"]])
            library.spectra_data["COLLISION_ENERGY"] = np.maximum(
                0, library.spectra_data["COLLISION_ENERGY"] + delta_ce
            )

        else:
            ce_alignment = alignment_library.spectra_data.groupby(by=["COLLISION_ENERGY"])["SPECTRAL_ANGLE"].mean()
            best_ce = ce_alignment.idxmax()
            pl.plot_mean_sa_ce(
                sa_ce_df=ce_alignment.to_frame().reset_index(),
                filename=results_dir / f"{spectra_file.stem}_mean_spectral_angle_ce.svg",
            )
            pl.plot_violin_sa_ce(
                sa_ce_df=alignment_library.spectra_data[["COLLISION_ENERGY", "SPECTRAL_ANGLE"]],
                filename=results_dir / f"{spectra_file.stem}_violin_spectral_angle_ce.svg",
            )
            library.spectra_data["COLLISION_ENERGY"] = best_ce
            with open(results_dir / f"{spectra_file.stem}_ce.txt", "w") as f:
                f.write(str(best_ce))
    else:
        best_ce = 35
        library.spectra_data["COLLISION_ENERGY"] = best_ce

        with open(results_dir / f"{spectra_file.stem}_ce.txt", "w") as f:
            f.write(str(best_ce))


def generate_spectral_lib(config_path: Union[str, Path]):
    """
    Create a SpectralLibrary object and generate the spectral library.

    # TODO full description
    :param config_path: path to config file
    :raises ValueError: spectral library output format is not supported as spectral library type
    """
    config = Config()
    config.read(config_path)
    library_input_type = config.library_input_type

    if library_input_type == "fasta":
        pp.digest(
            fasta=config.library_input,
            output=config.output,
            fragmentation=config.fragmentation,
            digestion=config.digestion,
            cleavages=config.cleavages,
            db=config.db,
            enzyme=config.enzyme,
            special_aas=config.special_aas,
            min_length=config.min_length,
            max_length=config.max_length,
        )
        library_file = config.output / "prosit_input.csv"
    elif library_input_type == "peptides":
        library_file = config.library_input
    else:
        raise ValueError(f'Library input type {library_input_type} not understood. Can only be "fasta" or "peptides".')
    spec_library = pp.gen_lib(library_file)

    pp_and_filter_step = ProcessStep(config.output, "speclib_filtered")

    if not pp_and_filter_step.is_done():
        spec_library = pp.process_and_filter_spectra_data(
            library=spec_library, model=config.models["intensity"], tmt_label=config.tag
        )
        spec_library.write_as_hdf5(config.output / "data" / f"{library_file.stem}_filtered.hdf5")
        pp_and_filter_step.mark_done()
    else:
        spec_library = Spectra.from_hdf5(config.output / "data" / f"{library_file.stem}_filtered.hdf5")

    no_of_spectra = len(spec_library.spectra_data)
    batchsize = 10000
    no_of_sections = np.math.ceil(no_of_spectra / batchsize)

    server_kwargs = {
        "server_url": config.prediction_server,
        "ssl": config.ssl,
    }

    spectral_library: Type[SpectralLibrary]
    results_path = config.output / "results"
    results_path.mkdir(exist_ok=True)

    if config.output_format == "msp":
        spectral_library = MSP
        out_file = results_path / "myPrositLib.msp"
    elif config.output_format == "spectronaut":
        spectral_library = Spectronaut
        out_file = results_path / "myPrositLib.csv"
    elif config.output_format == "dlib":
        spectral_library = DLib
        out_file = results_path / "myPrositLib.dlib"
    else:
        raise ValueError(f"{config.output_format} is not supported as spectral library type")

    if out_file.is_file():
        out_file.unlink()

    speclib = spectral_library(out_file, mode="w")

    spec_library.spectra_data.rename(
        columns={
            "MODIFIED_SEQUENCE": "peptide_sequences",
            "PRECURSOR_CHARGE": "precursor_charges",
            "COLLISION_ENERGY": "collision_energies",
            "FRAGMENTATION": "fragmentation_types",
        },
        inplace=True,
    )

    with Manager() as manager:

        # setup
        shared_queue = manager.Queue(maxsize=config.num_threads)
        prediction_progress = manager.Value("i", 0)
        writing_progress = manager.Value("i", 0)

        lock = manager.Lock()

        # Create a pool for producer processes
        pool = Pool(config.num_threads)

        try:
            for i in range(no_of_sections):
                pool.apply_async(
                    _make_predictions,
                    (
                        spec_library.spectra_data.iloc[i * batchsize : (i + 1) * batchsize],
                        config.models["intensity"],
                        config.models["irt"],
                        server_kwargs,
                        shared_queue,
                        prediction_progress,
                        lock,
                    ),
                )

            # Start the consumer process
            with tqdm(total=no_of_sections, desc="Getting predictions", disable=False) as predictor_pbar:
                with tqdm(total=no_of_sections, desc="Writing library", disable=False) as writer_pbar:
                    args = shared_queue, writing_progress
                    consumer_process = Process(target=speclib.async_write, args=args)
                    consumer_process.start()
                    while prediction_progress.value < no_of_sections:
                        time.sleep(0.5)
                        predictor_pbar.n = prediction_progress.value
                        predictor_pbar.refresh()
                        writer_pbar.n = writing_progress.value
                        writer_pbar.refresh()
                while writing_progress.value < no_of_sections:
                    time.sleep(0.5)
                    writer_pbar.n = writing_progress.value
                    writer_pbar.refresh()

            # Wait for all producer processes to finish
            pool.close()
            pool.join()

            # Signal the consumer process that producers have finished
            shared_queue.put(None)

            # Wait for the consumer process to finish
            consumer_process.join()

        except (KeyboardInterrupt, SystemExit):
            logger.error("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            pool.join()
            sys.exit(1)

        except Exception as e:
            logger.error("Caught Unknown exception, terminating workers")
            logger.error(traceback.format_exc())
            logger.error(e)
            pool.terminate()
            pool.join()
            sys.exit(1)

    logger.info("Finished writing the library to disk")


def _ce_calib(spectra_file: Path, config: Config) -> Spectra:
    ce_calib_step = ProcessStep(config.output, "ce_calib." + spectra_file.stem)
    if ce_calib_step.is_done():
        hdf5_path = config.output / "data" / spectra_file.with_suffix(".mzml.pred.hdf5").name
        if hdf5_path.is_file():
            library = Spectra.from_hdf5(hdf5_path)
            return library
        else:
            raise FileNotFoundError(f"{hdf5_path} not found but ce_calib.{spectra_file.stem} found. Please check.")
    library = _annotate_and_get_library(spectra_file, config)
    _get_best_ce(library, spectra_file, config)

    library.write_pred_as_hdf5(config.output / "data" / spectra_file.with_suffix(".mzml.pred.hdf5").name)

    ce_calib_step.mark_done()

    return library


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
    spectra_files = pp.list_spectra(input_dir=config.spectra, file_format=config.spectra_type)
    logger.info(f"Found {len(spectra_files)} files in the spectra directory.")

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


def _calculate_features(spectra_file: Path, config: Config):
    library = _ce_calib(spectra_file, config)

    calc_feature_step = ProcessStep(config.output, "calculate_features." + spectra_file.stem)
    if calc_feature_step.is_done():
        return

    server_kwargs = {
        "server_url": config.prediction_server,
        "ssl": config.ssl,
    }

    library.spectra_data.rename(
        columns={
            "MODIFIED_SEQUENCE": "peptide_sequences",
            "PRECURSOR_CHARGE": "precursor_charges",
            "COLLISION_ENERGY": "collision_energies",
            "FRAGMENTATION": "fragmentation_types",
        },
        inplace=True,
    )

    pred_intensities = pr.predict(
        library.spectra_data,
        model_name=config.models["intensity"],
        **server_kwargs,
    )

    pred_irts = pr.predict(library.spectra_data, model_name=config.models["irt"], **server_kwargs)

    library.spectra_data.rename(
        columns={
            "peptide_sequences": "MODIFIED_SEQUENCE",
            "precursor_charges": "PRECURSOR_CHARGE",
            "collision_energies": "COLLISION_ENERGY",
            "fragmentation_types": "FRAGMENTATION",
        },
        inplace=True,
    )

    library.add_matrix(pd.Series(pred_intensities["intensities"].tolist(), name="intensities"), FragmentType.PRED)
    library.add_column(pred_irts["irt"], name="PREDICTED_IRT")

    library.write_pred_as_hdf5(config.output / "data" / spectra_file.with_suffix(".mzml.pred.hdf5").name)

    # produce percolator tab files
    fdr_dir = config.output / "results" / config.fdr_estimation_method
    fdr_dir.mkdir(exist_ok=True)

    re.generate_features(
        library=library,
        search_type="original",
        output_file=fdr_dir / spectra_file.with_suffix(".original.tab").name,
        all_features=config.all_features,
        regression_method=config.curve_fitting_method,
    )
    re.generate_features(
        library=library,
        search_type="rescore",
        output_file=fdr_dir / spectra_file.with_suffix(".rescore.tab").name,
        all_features=config.all_features,
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


def run_rescoring(config_path: Union[str, Path]):
    """
    Create a ReScore object and run the rescoring.

    # TODO full description
    :param config_path: path to config file
    """
    config = Config()
    config.read(config_path)

    # load spectra file names
    spectra_files = pp.list_spectra(input_dir=config.spectra, file_format=config.spectra_type)
    logger.info(f"Found {len(spectra_files)} files in the spectra directory.")

    proc_dir = config.output / "proc"
    proc_dir.mkdir(parents=True, exist_ok=True)

    spectra_files = _preprocess(spectra_files, config)

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
