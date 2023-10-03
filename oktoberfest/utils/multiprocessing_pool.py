import logging
import signal
import sys
import traceback
import warnings
from multiprocessing import Pool, pool
from typing import List

from tqdm.auto import tqdm

logger = logging.getLogger(__name__)


class JobPool:
    """JobPool class for multiprocessing."""

    results: List[pool.AsyncResult]
    warning_filter: str
    pool: pool.Pool

    def __init__(self, processes: int = 1, warning_filter: str = "default"):
        """Initialize JobPool."""
        self.warning_filter = warning_filter
        self.pool = Pool(processes, self.init_worker)
        self.results = []

    def apply_async(self, f, args):
        """Apply async."""
        r = self.pool.apply_async(f, args)
        self.results.append(r)

    def init_worker(self):
        """Initialize the worker."""
        return init_worker(self.warning_filter)

    def check_pool(self):
        """Check the pool."""
        try:
            outputs = []
            res_len = len(self.results)
            with tqdm(total=res_len, desc="Waiting for tasks to complete", leave=True) as progress:
                for res in self.results:
                    outputs.append(res.get(timeout=10000))  # 10000 seconds = ~3 hours
                    progress.update(1)

            self.pool.close()
            self.pool.join()
            return outputs
        except (KeyboardInterrupt, SystemExit):
            logger.error("Caught KeyboardInterrupt, terminating workers")
            self.pool.terminate()
            self.pool.join()
            sys.exit(1)
        except Exception as e:
            logger.error("Caught Unknown exception, terminating workers")
            logger.error(traceback.format_exc())
            logger.error(e)
            self.pool.terminate()
            self.pool.join()
            sys.exit(1)


def init_worker(warning_filter):
    """Initialize worker given warning filter."""
    # set warning_filter for the child processes
    warnings.simplefilter(warning_filter)

    # causes child processes to ignore SIGINT signal and lets main process handle
    # interrupts instead (https://noswap.com/blog/python-multiprocessing-keyboardinterrupt)
    signal.signal(signal.SIGINT, signal.SIG_IGN)
