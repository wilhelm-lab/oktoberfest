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
                            logger.warning(
                                f"Task {label} timed out after {self.timeout/60:.1f} min — retrying ({retries+1}/{self.max_retries})."
                            )
                            next_round.append(self._resubmit(res))
                        else:
                            logger.error(f"Task {label} timed out after {self.max_retries} retries.")
                            status[label] = {"status": "timeout", "retries": retries}
                            progress.update(1)
                    except Exception as e:
                        retries = getattr(res, "_retries", 0)
                        if retries < self.max_retries:
                            logger.warning(
                                f"Task {label} failed with {e.__class__.__name__}: {e} — retrying ({retries+1}/{self.max_retries})."
                            )
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
        func = getattr(res, "_func")
        args = getattr(res, "_args")
        kwargs = getattr(res, "_kwargs")
        label = getattr(res, "_job_label", "unknown")
        tries = getattr(res, "_retries", 0) + 1

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
