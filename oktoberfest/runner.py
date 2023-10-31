import datetime
import json
import logging
from pathlib import Path
from typing import List, Type, Union

from spectrum_io.spectral_library import MSP, DLib, SpectralLibrary, Spectronaut

from oktoberfest import __copyright__, __version__
from oktoberfest import plotting as pl
from oktoberfest import predict as pr
from oktoberfest import preprocessing as pp
from oktoberfest import rescore as re

from .data.spectra import Spectra
from .utils import Config, JobPool, ProcessStep

logger = logging.getLogger(__name__)


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
        mzml_dir = config.output / "mzML"
        mzml_dir.mkdir(exist_ok=True)
        if spectra_file.suffix.lower() == ".raw":
            mzml_file = mzml_dir / spectra_file.with_suffix(".mzML").name
            pp.convert_spectra(spectra_file, mzml_file, thermo_exe=config.thermo_exe)
        else:
            mzml_file = spectra_file
        spectra = pp.load_spectra(mzml_file)
        search = pp.load_search(config.output / "msms" / spectra_file.with_suffix(".rescore").name)
        library = pp.merge_spectra_and_peptides(spectra, search)
        pp.annotate_spectral_library(library, mass_tol=config.mass_tolerance, unit_mass_tol=config.unit_mass_tolerance)
        library.write_as_hdf5(hdf5_path)  # write_metadata_annotation

    return library


def _get_best_ce(library: Spectra, spectra_file: Path, config: Config) -> int:
    results_dir = config.output / "results"
    results_dir.mkdir(exist_ok=True)
    if (library.spectra_data["FRAGMENTATION"] == "HCD").any():
        server_kwargs = {
            "url": config.prediction_server,
            "ssl": config.ssl,
            "intensity_model": config.models["intensity"],
            "irt_model": config.models["irt"],
        }
        alignment_library = pr.ce_calibration(library, **server_kwargs)
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
    else:
        best_ce = 35

    with open(results_dir / f"{spectra_file.stem}_ce.txt", "w") as f:
        f.write(str(best_ce))

    return best_ce


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
    spec_library = pp.process_and_filter_spectra_data(
        library=spec_library, model=config.models["intensity"], tmt_label=config.tag
    )

    no_of_spectra = len(spec_library.spectra_data)
    no_of_sections = no_of_spectra // 7000

    server_kwargs = {
        "url": config.prediction_server,
        "ssl": config.ssl,
        "intensity_model": config.models["intensity"],
        "irt_model": config.models["irt"],
        "job_type": "SpectralLibraryGeneration",
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

    for i in range(0, no_of_sections + 1):
        spectra_div = Spectra()
        if i < no_of_sections:
            spectra_div.spectra_data = spec_library.spectra_data.iloc[i * 7000 : (i + 1) * 7000]
            logger.info(f"Indices {i * 7000}, {(i + 1) * 7000}")
        elif (i * 7000) < no_of_spectra:
            spectra_div.spectra_data = spec_library.spectra_data.iloc[i * 7000 :]
            logger.info(f"Last Batch from index {i * 7000}")
            logger.info(f"Batch of size {len(spectra_div.spectra_data.index)}")
        else:
            break

        grpc_output_sec = pr.grpc_predict(spectra_div, **server_kwargs)

        out_lib = spectral_library(spectra_div.spectra_data, grpc_output_sec, out_file)
        out_lib.prepare_spectrum()
        out_lib.write()


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
    best_ce = _get_best_ce(library, spectra_file, config)

    library.spectra_data["COLLISION_ENERGY"] = best_ce
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

    processing_pool = JobPool(processes=config.num_threads)

    for spectra_file in spectra_files:
        processing_pool.apply_async(_ce_calib, [spectra_file, config])
    processing_pool.check_pool()


def _calculate_features(spectra_file: Path, config: Config):
    library = _ce_calib(spectra_file, config)

    calc_feature_step = ProcessStep(config.output, "calculate_features." + spectra_file.stem)
    if calc_feature_step.is_done():
        return

    server_kwargs = {
        "url": config.prediction_server,
        "ssl": config.ssl,
        "intensity_model": config.models["intensity"],
        "irt_model": config.models["irt"],
    }

    pr.grpc_predict(library, **server_kwargs)
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


def run_rescoring(config_path: Union[str, Path]):
    """
    Create a ReScore object and run the rescoring.

    # TODO full description
    :param config_path: path to config file
    :raises ValueError: if the provided fdr estimation method in the config is not recognized
    """
    config = Config()
    config.read(config_path)

    # load spectra file names
    spectra_files = pp.list_spectra(input_dir=config.spectra, file_format=config.spectra_type)
    logger.info(f"Found {len(spectra_files)} files in the spectra directory.")

    proc_dir = config.output / "proc"
    proc_dir.mkdir(parents=True, exist_ok=True)

    spectra_files = _preprocess(spectra_files, config)

    processing_pool = JobPool(processes=config.num_threads)

    for spectra_file in spectra_files:
        processing_pool.apply_async(_calculate_features, [spectra_file, config])
    processing_pool.check_pool()

    # rescoring

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
    job_type = conf.job_type

    # add file handler to root logger
    base_logger = logging.getLogger()
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(name)s::%(funcName)s %(message)s")
    suffix = datetime.datetime.now().strftime("%Y_%m_%d-%H_%M_%S")
    logging_output = conf.output / f"{job_type}_{suffix}.log"
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
