import logging
import os
from pathlib import Path

import spectrum_fundamentals.constants as c
from spectrum_fundamentals.fragments import compute_peptide_mass
from spectrum_fundamentals.mod_string import internal_without_mods, maxquant_to_internal
from spectrum_io import Spectronaut
from spectrum_io.spectral_library import MSP

from .ce_calibration import CeCalibration, SpectralLibrary
from .data.spectra import Spectra
from .re_score import ReScore
from .utils.config import Config

__version__ = "0.1.0"
__copyright__ = """Copyright (c) 2020-2021 Oktoberfest dev-team. All rights reserved.
Written by
- Wassim Gabriel (wassim.gabriel@tum.de),
- Ludwig Lautenbacher (ludwig.lautenbacher@tum.de),
- Matthew The (matthew.the@tum.de),
- Mario Picciani (mario.picciani@in.tum.de),
- Firas Hamood (firas.hamood@tum.de),
- Cecilia Jensen (cecilia.jensen@tum.de)
at the Technical University of Munich."""

logger = logging.getLogger(__name__)


def generate_spectral_lib(search_dir: str, config_path: str):
    """
    Create a SpectralLibrary object and generate the spectral library.

    :param search_dir: path to directory containing the msms.txt and raw files
    :param config_path: path to config file
    :raises ValueError: spectral library output format is not supported as spectral library type
    """
    spec_library = SpectralLibrary(path=search_dir, out_path=search_dir, config_path=config_path)
    spec_library.gen_lib()
    spec_library.library.spectra_data["MODIFIED_SEQUENCE"] = spec_library.library.spectra_data[
        "MODIFIED_SEQUENCE"
    ].apply(lambda x: "_" + x + "_")
    models_dict = spec_library.config.models
    spec_library.library.spectra_data["MODIFIED_SEQUENCE"] = maxquant_to_internal(
        spec_library.library.spectra_data["MODIFIED_SEQUENCE"], fixed_mods={}
    )
    spec_library.library.spectra_data["SEQUENCE"] = internal_without_mods(
        spec_library.library.spectra_data["MODIFIED_SEQUENCE"]
    )
    spec_library.library.spectra_data["PEPTIDE_LENGTH"] = spec_library.library.spectra_data["SEQUENCE"].apply(
        lambda x: len(x)
    )

    logger.info(f"No of sequences before Filtering is {len(spec_library.library.spectra_data['PEPTIDE_LENGTH'])}")
    spec_library.library.spectra_data = spec_library.library.spectra_data[
        (spec_library.library.spectra_data["PEPTIDE_LENGTH"] <= 30)
    ]
    spec_library.library.spectra_data = spec_library.library.spectra_data[
        (~spec_library.library.spectra_data["MODIFIED_SEQUENCE"].str.contains(r"\(ac\)"))
    ]
    spec_library.library.spectra_data = spec_library.library.spectra_data[
        (~spec_library.library.spectra_data["MODIFIED_SEQUENCE"].str.contains(r"\(Acetyl \(Protein N-term\)\)"))
    ]
    spec_library.library.spectra_data = spec_library.library.spectra_data[
        (~spec_library.library.spectra_data["SEQUENCE"].str.contains("U"))
    ]
    spec_library.library.spectra_data = spec_library.library.spectra_data[
        spec_library.library.spectra_data["PRECURSOR_CHARGE"] <= 6
    ]
    spec_library.library.spectra_data = spec_library.library.spectra_data[
        spec_library.library.spectra_data["PEPTIDE_LENGTH"] >= 7
    ]
    logger.info(f"No of sequences after Filtering is {len(spec_library.library.spectra_data['PEPTIDE_LENGTH'])}")

    tmt_model = False
    for _, value in models_dict.items():
        if value:
            if "TMT" in value:
                tmt_model = True
    if tmt_model and spec_library.config.tag != "":
        unimod_tag = c.TMT_MODS[spec_library.config.tag]
        spec_library.library.spectra_data["MODIFIED_SEQUENCE"] = maxquant_to_internal(
            spec_library.library.spectra_data["MODIFIED_SEQUENCE"],
            fixed_mods={"C": "C[UNIMOD:4]", "^_": f"_{unimod_tag}", "K": f"K{unimod_tag}"},
        )
    else:
        spec_library.library.spectra_data["MODIFIED_SEQUENCE"] = maxquant_to_internal(
            spec_library.library.spectra_data["MODIFIED_SEQUENCE"]
        )
    spec_library.library.spectra_data["MASS"] = spec_library.library.spectra_data["MODIFIED_SEQUENCE"].apply(
        lambda x: compute_peptide_mass(x)
    )
    no_of_spectra = len(spec_library.library.spectra_data)
    no_of_sections = no_of_spectra // 7000
    for i in range(0, no_of_sections + 1):
        spectra_div = Spectra()
        if i < no_of_sections:
            spectra_div.spectra_data = spec_library.library.spectra_data.iloc[i * 7000 : (i + 1) * 7000]
            logger.info(f"Indices {i * 7000}, {(i + 1) * 7000}")
        elif (i * 7000) < no_of_spectra:
            spectra_div.spectra_data = spec_library.library.spectra_data.iloc[i * 7000 :]
            logger.info(f"Last Batch from index {i * 7000}")
            logger.info(f"Batch of size {len(spectra_div.spectra_data.index)}")
        else:
            break

        grpc_output_sec = spec_library.grpc_predict(spectra_div)
        if spec_library.config.output_format == "msp":
            out_lib_msp = MSP(
                spectra_div.spectra_data, grpc_output_sec, os.path.join(spec_library.results_path, "myPrositLib.msp")
            )
            out_lib_msp.prepare_spectrum()
            out_lib_msp.write()
        elif spec_library.config.output_format == "spectronaut":
            out_lib_spectronaut = Spectronaut(
                spectra_div.spectra_data, grpc_output_sec, os.path.join(spec_library.results_path, "myPrositLib.csv")
            )
            out_lib_spectronaut.prepare_spectrum()
            out_lib_spectronaut.write()
        else:
            raise ValueError(f"{spec_library.config.output_format} is not supported as spectral library type")


def run_ce_calibration(msms_path: str, search_dir: str, config_path: str):
    """
    Create a CeCalibration object and run the CE calibration.

    :param msms_path: path to msms folder
    :param search_dir: path to directory containing the msms.txt and raw files
    :param config_path: path to config file
    :raises ValueError: raw_type is not supported as rawfile-type
    """
    ce_calib = CeCalibration(search_path=msms_path, raw_path=search_dir, out_path=search_dir, config_path=config_path)
    df_search = ce_calib._load_search()
    raw_type = ce_calib.config.raw_type
    if raw_type == "thermo":
        extension = ".raw"
    elif raw_type == "mzml":
        extension = ".mzml"
    else:
        raise ValueError(f"{raw_type} is not supported as rawfile-type")

    ce_calib.raw_path = Path(
        os.path.join(
            ce_calib.raw_path,
            [os.path.basename(f) for f in os.listdir(ce_calib.raw_path) if f.lower().endswith(extension)][0],
        )
    )
    ce_calib.perform_alignment(df_search)
    with open(os.path.join(ce_calib.results_path, "ce.txt"), "w") as f:
        f.write(str(ce_calib.best_ce))


def run_rescoring(msms_path: str, search_dir: str, config_path: str):
    """
    Create a ReScore object and run the rescoring.

    :param msms_path: path to msms folder
    :param search_dir: path to directory containing the msms.txt and raw files
    :param config_path: path to config file
    """
    re_score = ReScore(search_path=msms_path, raw_path=search_dir, out_path=search_dir, config_path=config_path)
    re_score.get_raw_files()
    re_score.split_msms()
    re_score.calculate_features()

    re_score.merge_input("rescore")
    re_score.merge_input("original")

    re_score.rescore_with_perc("rescore")
    re_score.rescore_with_perc("original")


def run_job(search_dir: str, config_path: str):
    """
    Run oktoberfest based on job type given in the config file.

    :param search_dir: path to directory containing the msms.txt and raw files
    :param config_path: path to config file as a string
    :raises ValueError: In case the job_type in the provided config file is not known
    """
    if not config_path:
        config_path = os.path.join(search_dir, "config.json")
    conf = Config()
    conf.read(config_path)
    job_type = conf.job_type
    if conf.search_path:
        msms_path = conf.search_path
    else:
        msms_path = os.path.join(search_dir, "msms.txt")

    if job_type == "SpectralLibraryGeneration":
        generate_spectral_lib(search_dir, config_path)
    elif job_type == "CollisionEnergyCalibration":
        run_ce_calibration(msms_path, search_dir, config_path)
    elif job_type == "Rescoring":
        run_rescoring(msms_path, search_dir, config_path)
    else:
        raise ValueError(f"Unknown job_type in config: {job_type}")
