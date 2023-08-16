from pathlib import Path

from oktoberfest import runner


def integration_test():
    """Integration test for the rescoring."""
    runner.run_job(Path("data/plasma/config.json"))


if __name__ == "__main__":
    integration_test()
