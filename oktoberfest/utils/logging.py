import contextlib
import os
import warnings


@contextlib.contextmanager
def mute_stdout(ignore_warnings: bool = False):
    """Mute print statements and user warnings from packages that aren't properly handled."""
    with open(os.devnull, "w") as devnull, contextlib.redirect_stdout(devnull):
        if ignore_warnings:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                yield
        else:
            yield
