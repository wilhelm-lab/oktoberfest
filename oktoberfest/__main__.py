import argparse

from rich import traceback

from oktoberfest import __version__, runner

"""triqler.__main__: executed when bootstrap directory is called as script."""

"""Command-line interface."""
# import click


def _parse_args():
    """Parse search_dir and config_path arguments."""
    apars = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)  # , add_help=False)

    # option_group = apars.add_argument_group("optional arguments")
    # option_group.add_argument("-h", "--help", action="help", help="show this help message and exit")

    apars.add_argument("-v", "--version", action="version", version=f"{__version__}")
    apars.add_argument(
        "-c",
        "--config_path",
        default=None,
        metavar="CONFIG",
        required=True,
        help=(
            "Path to config file in json format."
            "If this argument is not specified, we try to find and use a file called config.json in <search_dir>."
        ),
    )

    args = apars.parse_args()
    return args


def main():
    """Execution of oktoberfest from terminal."""
    args = _parse_args()
    runner.run_job(args.config_path)


if __name__ == "__main__":
    traceback.install()
    main()  # pragma: no cover
