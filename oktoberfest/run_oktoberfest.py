import argparse
import os
import sys

from oktoberfest import __copyright__, __version__, logger, runner


def parse_args():
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
    logger.info(f"Oktoberfest version {__version__}\n{__copyright__}")
    logger.info(f'Issued command: {os.path.basename(__file__)} {" ".join(map(str, sys.argv[1:]))}')

    args = parse_args()
    runner.run_job(args.config_path)


if __name__ == "__main__":
    main()
