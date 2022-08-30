"""triqler.__main__: executed when bootstrap directory is called as script."""

"""Command-line interface."""
import click
from rich import traceback

from .oktoberfest import main

# main()

if __name__ == "__main__":
    traceback.install()
    main(prog_name="batchglm")  # pragma: no cover
