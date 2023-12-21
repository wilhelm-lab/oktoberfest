import argparse

from rich import traceback

from oktoberfest import runner

"""triqler.__main__: executed when bootstrap directory is called as script."""

"""Command-line interface."""
# import click


def _parse_args():
    """Parse search_dir and config_path arguments."""
    apars = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    apars.add_argument(
        "--config_path",
        default=None,
        metavar="C",
        help="""Path to config file in json format. \\
                If this argument is not specified, we try to find and use a file called config.json in <search_dir>.""",
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
