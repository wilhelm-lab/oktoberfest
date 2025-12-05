import logging
import signal
import sys
import traceback
import warnings
from multiprocessing import Pool, TimeoutError

from tqdm.auto import tqdm

logger = logging.getLogger(__name__)


# class JobPool:
#     """JobPool class for multiprocessing."""

#     results: list[pool.AsyncResult]
#     warning_filter: str
#     pool: pool.Pool

#     def __init__(self, processes: int = 1, warning_filter: str = "default"):
#         """Initialize JobPool."""
#         self.warning_filter = warning_filter
#         self.pool = Pool(processes, self.init_worker)
#         self.results = []

#     def apply_async(self, f, args, **kwargs):
#         """Apply async."""
#         r = self.pool.apply_async(f, args=args, kwds=kwargs)
#         self.results.append(r)

#     def init_worker(self):
#         """Initialize the worker."""
#         return init_worker(self.warning_filter)

#     def check_pool(self):
#         """Check the pool."""
#         try:
#             outputs = []
#             res_len = len(self.results)
#             with tqdm(total=res_len, desc="Waiting for tasks to complete", leave=True) as progress:
#                 for res in self.results:
#                     outputs.append(res.get(timeout=10000))  # 10000 seconds = ~3 hours
#                     progress.update(1)
#             self.pool.close()
#             self.pool.join()
#             return outputs
#         except (KeyboardInterrupt, SystemExit):
#             logger.error("Caught KeyboardInterrupt, terminating workers")
#             self.pool.terminate()
#             self.pool.join()
#             sys.exit(1)
#         except Exception as e:
#             logger.error("Caught Unknown exception, terminating workers")
#             logger.error(traceback.format_exc())
#             logger.error(e)
#             self.pool.terminate()
#             self.pool.join()
#             sys.exit(1)


# def init_worker(warning_filter):
#     """Initialize worker given warning filter."""
#     # set warning_filter for the child processes
#     warnings.simplefilter(warning_filter)

#     # causes child processes to ignore SIGINT signal and lets main process handle
#     # interrupts instead (https://noswap.com/blog/python-multiprocessing-keyboardinterrupt)

#     signal.signal(signal.SIGINT, signal.SIG_IGN)
class JobPool:
    """
    Multiprocessing pool for side-effect-only workers (no return values).
    Tracks per-job status instead of collecting results.
    """

    def __init__(self, processes: int = 1, warning_filter: str = "default",
                 timeout: int = 1200, max_retries: int = 2):
        self.warning_filter = warning_filter
        self.timeout = timeout          # 20 minutes
        self.max_retries = max_retries
        self.pool = Pool(processes, self._init_worker)
        self.results = []

    def _init_worker(self):
        _init_worker(self.warning_filter)

    def apply_async(self, func, args=(), **kwargs):
        """
        Submit a job. Workers may not return anything; we only track status.
        You can pass job_label=... in kwargs for clearer logs.
        """
        job_label = kwargs.pop("job_label", args[0] if args else f"job_{len(self.results)}")
        r = self.pool.apply_async(func, args=args, kwds=kwargs)
        # attach metadata for tracking
        r._job_label = str(job_label)
        r._func = func
        r._args = args
        r._kwargs = kwargs
        r._retries = 0
        r._index = len(self.results)
        self.results.append(r)
        return r

    def check_pool(self):
        """
        Waits for all jobs, retrying failures/timeouts up to max_retries.
        Returns a dict: {job_label: {"status": "ok"|"timeout"|"error", "retries": n}}
        """
        status = {}
        remaining = list(self.results)

        with tqdm(total=len(remaining), desc="Processing tasks", leave=True) as progress:
            while remaining:
                next_round = []
                for res in remaining:
                    label = getattr(res, "_job_label", "unknown")
                    try:
                        # We ignore the return value on purpose (workers have side effects)
                        res.get(timeout=self.timeout)
                        status[label] = {"status": "ok", "retries": getattr(res, "_retries", 0)}
                        logger.debug(f"Task {label} completed (no return).")
                        progress.update(1)
                    except TimeoutError:
                        retries = getattr(res, "_retries", 0)
                        if retries < self.max_retries:
                            logger.warning(f"Task {label} timed out after {self.timeout/60:.1f} min — retrying ({retries+1}/{self.max_retries}).")
                            next_round.append(self._resubmit(res))
                        else:
                            logger.error(f"Task {label} timed out after {self.max_retries} retries.")
                            status[label] = {"status": "timeout", "retries": retries}
                            progress.update(1)
                    except Exception as e:
                        retries = getattr(res, "_retries", 0)
                        if retries < self.max_retries:
                            logger.warning(f"Task {label} failed with {e.__class__.__name__}: {e} — retrying ({retries+1}/{self.max_retries}).")
                            logger.debug(traceback.format_exc())
                            next_round.append(self._resubmit(res))
                        else:
                            logger.error(f"Task {label} failed after {self.max_retries} retries.")
                            logger.error(traceback.format_exc())
                            status[label] = {"status": "error", "retries": retries}
                            progress.update(1)
                remaining = next_round

        self.pool.close()
        self.pool.join()
        logger.info("All tasks processed.")
        return status

    def _resubmit(self, res):
        func   = getattr(res, "_func")
        args   = getattr(res, "_args")
        kwargs = getattr(res, "_kwargs")
        label  = getattr(res, "_job_label", "unknown")
        tries  = getattr(res, "_retries", 0) + 1

        new_res = self.pool.apply_async(func, args=args, kwds=kwargs)
        new_res._job_label = label
        new_res._func = func
        new_res._args = args
        new_res._kwargs = kwargs
        new_res._retries = tries
        new_res._index = getattr(res, "_index", 0)
        logger.debug(f"Resubmitted task {label} (retry {tries}).")
        return new_res

def _init_worker(warning_filter):
    warnings.simplefilter(warning_filter)
    signal.signal(signal.SIGINT, signal.SIG_IGN)