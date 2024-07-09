import datetime
import json
import logging
from pathlib import Path
from typing import Union

from oktoberfest import __copyright__, __version__
from oktoberfest.jobs import ce_calibration, rescoring, spectral_library

from .utils import Config

logger = logging.getLogger(__name__)


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
            spectral_library.generate_spectral_lib(config_path)
        elif job_type == "CollisionEnergyCalibration":
            ce_calibration.run_ce_calibration(config_path)
        elif job_type == "Rescoring":
            rescoring.run_rescoring(config_path)
        else:
            raise ValueError(f"Unknown job_type in config: {job_type}")
    finally:
        file_handler.close()
        base_logger.removeHandler(file_handler)
