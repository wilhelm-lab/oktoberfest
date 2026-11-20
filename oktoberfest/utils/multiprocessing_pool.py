import logging
import multiprocessing
import signal
import sys
import traceback
import warnings
from concurrent.futures import Future, ProcessPoolExecutor, as_completed

from tqdm.auto import tqdm

logger = logging.getLogger(__name__)


def mp_context() -> multiprocessing.context.BaseContext:
    """
    Get the start method of Oktoberfest's worker processes: forkserver where available (Linux, macOS), else spawn.

    A forked worker copies the parent with all its locks, also those held by other threads at that moment (e.g. of
    the progress bars or of the Koina client), which can deadlock the worker (Python warns about this since 3.12 and
    no longer forks by default from 3.14 on). Forkserver and spawn start every worker from a fresh, single-threaded
    process. Scripts that start Oktoberfest jobs therefore need an ``if __name__ == "__main__":`` guard.

    :return: the multiprocessing context
    """
    methods = multiprocessing.get_all_start_methods()
    return multiprocessing.get_context("forkserver" if "forkserver" in methods else "spawn")


def log_files() -> list[tuple[str, str]]:
    """Get the file and format of every file handler of the root logger, so that workers can log to the same files."""
    return [
        (handler.baseFilename, handler.formatter._fmt if handler.formatter else "%(message)s")
        for handler in logging.getLogger().handlers
        if isinstance(handler, logging.FileHandler)
    ]


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
            mp_context=mp_context(),
            initializer=init_worker,
            initargs=(self.warning_filter, log_files()),
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
                    outputs[future_to_idx[future]] = future.result(timeout=100000)  # ~28 hours
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


def init_worker(warning_filter: str = "default", files: list[tuple[str, str]] | None = None):
    """
    Initialize a worker process: warning filter, log files of the parent and no SIGINT handling.

    :param warning_filter: warning filter of the child processes
    :param files: (file name, format) of the parent's log files (see :func:`log_files`); workers that are not forked
        do not inherit the parent's handlers, so their messages are appended to the same files
    """
    # set warning_filter for the child processes
    warnings.simplefilter(warning_filter)

    root = logging.getLogger()
    known = {handler.baseFilename for handler in root.handlers if isinstance(handler, logging.FileHandler)}
    for file_name, log_format in files or []:
        if file_name not in known:
            handler = logging.FileHandler(file_name)
            handler.setLevel(logging.DEBUG)
            handler.setFormatter(logging.Formatter(log_format))
            root.addHandler(handler)

    # causes child processes to ignore SIGINT signal and lets main process handle
    # interrupts instead (https://noswap.com/blog/python-multiprocessing-keyboardinterrupt)
    signal.signal(signal.SIGINT, signal.SIG_IGN)
