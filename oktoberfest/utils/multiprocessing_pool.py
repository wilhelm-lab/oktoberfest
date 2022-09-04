import logging
import signal
import sys
import traceback
import warnings
from multiprocessing import Pool

logger = logging.getLogger(__name__)


class JobPool:
    """JobPool class for multiprocessing."""

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

    def check_pool(self, print_progress_every: int = -1):
        """Check the pool."""
        try:
            outputs = list()
            for res in self.results:
                outputs.append(res.get(timeout=10000))  # 10000 seconds = ~3 hours
                if print_progress_every > 0 and len(outputs) % print_progress_every == 0:
                    logger.info(
                        f' {len(outputs)} / {len(self.results)} {"%.2f" % (float(len(outputs)) / len(self.results) * 100)}%'
                    )
            self.pool.close()
            self.pool.join()
            return outputs
        except (KeyboardInterrupt, SystemExit):
            logger.error("Caught KeyboardInterrupt, terminating workers")
            self.pool.terminate()
            self.pool.join()
            sys.exit()
        except Exception as e:
            logger.error("Caught Unknown exception, terminating workers")
            logger.error(traceback.print_exc())
            logger.error(e)
            self.pool.terminate()
            self.pool.join()
            sys.exit()


def init_worker(warning_filter):
    """Initialize worker given warning filter."""
    # set warning_filter for the child processes
    warnings.simplefilter(warning_filter)

    # causes child processes to ignore SIGINT signal and lets main process handle
    # interrupts instead (https://noswap.com/blog/python-multiprocessing-keyboardinterrupt)
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def add_one(i: int) -> int:
    """Add 1 to i."""
    return i + 1


def unit_test():
    """Unit test for multiprocessing."""
    pool = JobPool(4)
    for i in range(20):
        pool.apply_async(add_one, [i])
    results = pool.check_pool()
    print(results)
