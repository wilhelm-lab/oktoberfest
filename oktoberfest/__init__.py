"""Initialize logger."""
__version__ = "0.1.0"
import logging.handlers
import sys
import time

CONSOLE_LOG_LEVEL = logging.INFO
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
if len(logger.handlers) == 0:
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(name)s::%(funcName)s %(message)s")
    formatter.converter = time.gmtime
    # add console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(CONSOLE_LOG_LEVEL)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # add error handler
    error_handler = logging.StreamHandler()
    error_handler.setLevel(logging.ERROR)
    error_handler.setFormatter(formatter)
    logger.addHandler(error_handler)
else:
    logger.info("Logger already initizalized. Resuming normal operation.")
