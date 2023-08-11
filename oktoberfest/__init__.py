"""Initialize logger."""

__version__ = "0.4.0"
__copyright__ = """Copyright (c) 2020-2021 Oktoberfest dev-team. All rights reserved.
Written by
- Wassim Gabriel (wassim.gabriel@tum.de),
- Ludwig Lautenbacher (ludwig.lautenbacher@tum.de),
- Matthew The (matthew.the@tum.de),
- Mario Picciani (mario.picciani@in.tum.de),
- Firas Hamood (firas.hamood@tum.de),
- Cecilia Jensen (cecilia.jensen@tum.de)
at the Technical University of Munich."""
import logging.handlers
import sys
import time

from . import runner

CONSOLE_LOG_LEVEL = logging.INFO
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
if len(logger.handlers) == 0:
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(name)s::%(funcName)s %(message)s")
    converter = time.gmtime
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
