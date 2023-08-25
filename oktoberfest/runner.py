import logging
import os
from pathlib import Path
from typing import Union

from spectrum_io import Spectronaut
from spectrum_io.spectral_library import MSP

from oktoberfest import plotting as pl
from oktoberfest import predict as pr
from oktoberfest import preprocessing as pp
from oktoberfest import rescore as re

from .data.spectra import Spectra
from .utils.config import Config
from .utils.process_step import ProcessStep

__version__ = "0.1.0"
__copyright__ = """Copyright (c) 2020-2021 Oktoberfest dev-team. All rights reserved.
Written by
- Wassim Gabriel (wassim.gabriel@tum.de),
- Ludwig Lautenbacher (ludwig.lautenbacher@tum.de),
- Matthew The (matthew.the@tum.de),
- Mario Picciani (mario.picciani@tum.de),
- Firas Hamood (firas.hamood@tum.de),
- Cecilia Jensen (cecilia.jensen@tum.de)
at the Technical University of Munich."""

logger = logging.getLogger(__name__)


def generate_spectral_lib(search_dir: Union[str, Path], config_path: Union[str, Path]):
    """
    Create a SpectralLibrary object and generate the spectral library.

    :param search_dir: path to directory containing the msms.txt and raw files
    :param config_path: path to config file
    :raises ValueError: spectral library output format is not supported as spectral library type
    """
    config = Config()
    config.read(config_path)
    spec_library = Spectra()
    spec_library = pp.gen_lib(config, spec_library)
    spec_library_filtered = pp.process_and_filter_spectra_data(config, spec_library)

    no_of_spectra = len(spec_library_filtered.spectra_data)
    no_of_sections = no_of_spectra // 7000
    for i in range(0, no_of_sections + 1):
        spectra_div = Spectra()
        if i < no_of_sections:
            spectra_div.spectra_data = spec_library_filtered.spectra_data.iloc[i * 7000 : (i + 1) * 7000]
            logger.info(f"Indices {i * 7000}, {(i + 1) * 7000}")
        elif (i * 7000) < no_of_spectra:
            spectra_div.spectra_data = spec_library_filtered.spectra_data.iloc[i * 7000 :]
            logger.info(f"Last Batch from index {i * 7000}")
            logger.info(f"Batch of size {len(spectra_div.spectra_data.index)}")
        else:
            break

        grpc_output_sec = pr.grpc_predict(config, spectra_div)
        results_path = config.output / "results"
        results_path.mkdir(exist_ok=True)
        if config.output_format == "msp":
            out_lib_msp = MSP(spectra_div.spectra_data, grpc_output_sec, os.path.join(results_path, "myPrositLib.msp"))
            out_lib_msp.prepare_spectrum()
            out_lib_msp.write()
        elif config.output_format == "spectronaut":
            out_lib_spectronaut = Spectronaut(
                spectra_div.spectra_data, grpc_output_sec, os.path.join(results_path, "myPrositLib.csv")
            )
            out_lib_spectronaut.prepare_spectrum()
            out_lib_spectronaut.write()
        else:
            raise ValueError(f"{config.output_format} is not supported as spectral library type")


def run_ce_calibration(
    msms_path: Union[str, Path],
    search_dir: Union[str, Path],
    config_path: Union[str, Path],
    glob_pattern: str,
    output_path: Union[str, Path],
):
    """
    Create a CeCalibration object and run the CE calibration.

    :param msms_path: path to msms folder
    :param search_dir: path to directory containing the msms.txt and raw files
    :param config_path: path to config file
    :param glob_pattern: the pattern for raw file extensions to search for in search_dir
    :param output_path: path to the output folder if specified in the config file, else search_dir
    """
    if isinstance(search_dir, str):
        search_dir = Path(search_dir)

    config = Config()
    config.read(config_path)
    library = Spectra()

    spectra_type = config.spectra_type
    if spectra_type == "raw":
        glob_pattern = "*.[rR][aA][wW]"
    elif spectra_type == "mzml":
        glob_pattern = "*.[mM][zZ][mM][lL]"

    for raw_file in search_dir.glob(glob_pattern):
        df = pp.load_search(config)
        best_ce = pr.perform_alignment(config, library, df)
        with open(config.output / "results" / f"{raw_file.stem}_ce.txt", "w") as f:
            f.write(str(best_ce))


def run_rescoring(
    msms_path: Union[str, Path],
    search_dir: Union[str, Path],
    config_path: Union[str, Path],
    output_path: Union[str, Path],
):
    """
    Create a ReScore object and run the rescoring.

    :param msms_path: path to msms folder
    :param search_dir: path to directory containing the msms.txt and raw files
    :param config_path: path to config file
    :param output_path: path to the output folder if specified in the config file, else search_dir
    """
    logger.info("Starting rescoring run...")
    # re_score = ReScore(search_path=msms_path, raw_path=search_dir, out_path=output_path, config_path=config_path)
    config = Config()
    config.read(config_path)

    # re_score.get_raw_files()
    raw_files = pp.get_raw_files(config)
    # re_score.split_msms()
    pp.split_msms(config, ProcessStep(config.output, "split_msms"))
    # re_score.calculate_features()
    re.calculate_features(config, Path(config_path), raw_files)

    # re_score.merge_input("rescore")
    re.merge_input(
        config,
        ProcessStep(config.output, "merge_input_prosit"),
        ProcessStep(config.output, "merge_input_andromeda"),
        raw_files,
        "rescore",
    )
    # re_score.merge_input("original")
    re.merge_input(
        config,
        ProcessStep(config.output, "merge_input_prosit"),
        ProcessStep(config.output, "merge_input_andromeda"),
        raw_files,
        "original",
    )

    # re_score.rescore("rescore")
    re.rescore(
        config,
        ProcessStep(config.output, "percolator_prosit"),
        ProcessStep(config.output, "percolator_andromeda"),
        "rescore",
    )
    # re_score.rescore("original")
    re.rescore(
        config,
        ProcessStep(config.output, "percolator_prosit"),
        ProcessStep(config.output, "percolator_andromeda"),
        "original",
    )
    # plot_all(re_score.get_percolator_folder_path(), re_score.config.fdr_estimation_method)
    fdr_output = config.output / "results" / config.fdr_estimation_method
    pl.plot_all(fdr_output, config.fdr_estimation_method)

    logger.info("Finished rescoring.")


def run_job(config_path: Union[str, Path]):
    """
    Run oktoberfest based on job type given in the config file.

    :param config_path: path to config file as a string
    :raises ValueError: In case the job_type in the provided config file is not known
    """
    if not config_path:
        config_path = "./config.json"
    if isinstance(config_path, str):
        config_path = Path(config_path)
    conf = Config()
    conf.read(config_path)
    job_type = conf.job_type
    search_dir = conf.spectra
    output_path = conf.output
    msms_path = conf.search_results
    if job_type == "SpectralLibraryGeneration":
        generate_spectral_lib(search_dir, config_path)
    elif job_type == "CollisionEnergyCalibration":
        spectra_type = conf.spectra_type
        if spectra_type == "raw":
            glob_pattern = "*.[rR][aA][wW]"
        elif spectra_type == "mzml":
            glob_pattern = "*.[mM][zZ][mM][lL]"
        else:
            raise ValueError(f"{spectra_type} is not supported as spectra type (supported: raw and mzml).")
        run_ce_calibration(msms_path, search_dir, config_path, glob_pattern, output_path)
    elif job_type == "Rescoring":
        run_rescoring(msms_path, search_dir, config_path, output_path)
    else:
        raise ValueError(f"Unknown job_type in config: {job_type}")
