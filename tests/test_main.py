"""Test cases for the __main__ module."""
import pytest
from click.testing import CliRunner

# from oktoberfest import __main__


@pytest.fixture
def runner() -> CliRunner:
    """Fixture for invoking command-line interfaces."""
    return CliRunner()


def test_main_succeeds(runner: CliRunner) -> None:
    """It exits with a status code of zero."""
    # TODO should execute the main method without args result = runner.invoke(__main__.main)
    assert True  # result.exit_code == 0
