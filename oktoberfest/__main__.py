from rich import traceback

from .oktoberfest import main

"""triqler.__main__: executed when bootstrap directory is called as script."""

"""Command-line interface."""
# import click

# main()

if __name__ == "__main__":
    traceback.install()
    main(prog_name="batchglm")  # pragma: no cover
