import logging
import signal
import sys
import traceback
import warnings
from concurrent.futures import Future, ProcessPoolExecutor, as_completed

from tqdm.auto import tqdm

logger = logging.getLogger(__name__)


class JobPool:
    """Small wrapper around ProcessPoolExecutor."""

    results: list[Future]
    warning_filter: str
    executor: ProcessPoolExecutor

    def __init__(self, processes: int = 1, warning_filter: str = "default"):
        """Initialize JobPool."""
        self.warning_filter = warning_filter
        self.executor = ProcessPoolExecutor(
            max_workers=processes,
            initializer=init_worker,
            initargs=(self.warning_filter,),
        )
        self.results = []

    def apply_async(self, f, args, **kwargs):
        """Apply async."""
        future = self.executor.submit(f, *args, **kwargs)
        self.results.append(future)
        return future

    def check_pool(self):
        """Check the pool."""
        try:
            outputs = [None] * len(self.results)
            future_to_idx = {future: idx for idx, future in enumerate(self.results)}
            res_len = len(self.results)
            with tqdm(total=res_len, desc="Waiting for tasks to complete", leave=True) as progress:
                for future in as_completed(self.results):
                    outputs[future_to_idx[future]] = future.result(timeout=10000)  # ~3 hours
                    progress.update(1)
            self.executor.shutdown(wait=True, cancel_futures=False)
            return outputs
        except (KeyboardInterrupt, SystemExit):
            logger.error("Caught KeyboardInterrupt, terminating workers")
            self.executor.shutdown(wait=False, cancel_futures=True)
            sys.exit(1)
        except Exception as e:
            logger.error("Caught Unknown exception, terminating workers")
            logger.error(traceback.format_exc())
            logger.error(e)
            self.executor.shutdown(wait=False, cancel_futures=True)
            sys.exit(1)


def init_worker(warning_filter):
    """Initialize worker given warning filter."""
    # set warning_filter for the child processes
    warnings.simplefilter(warning_filter)

    # causes child processes to ignore SIGINT signal and lets main process handle
    # interrupts instead (https://noswap.com/blog/python-multiprocessing-keyboardinterrupt)
    signal.signal(signal.SIGINT, signal.SIG_IGN)
