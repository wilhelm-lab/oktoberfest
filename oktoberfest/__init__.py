import logging.handlers
import time
from .ce_calibration import CeCalibration
from .re_score import ReScore


CONSOLE_LOG_LEVEL = logging.INFO
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
if len(logger.handlers) == 0:
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(name)s::%(funcName)s %(message)s")
    formatter.converter = time.gmtime
    # add console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(CONSOLE_LOG_LEVEL)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
else:
    logger.info('Logger already initizalized. Resuming normal operation.')
