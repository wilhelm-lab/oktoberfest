import shutil
from pathlib import Path

import pandas as pd
import pytest

from oktoberfest import pp
from oktoberfest.utils import Config


class TestCustomMods:
    """Test class for testing feature of custom modifications."""

    def test_mods_maxquant(self, exp_data):
        """Test of custom modifications in Maxquant results."""
        config = Config()
        config.read(Path(__file__).parent / "configs/mods.json")

        msms_output = config.output / "msms"
        msms_output.mkdir(exist_ok=True, parents=True)
        internal_search_file = msms_output / "msms.prosit"
        tmt_label = config.tag
        search_results = pp.convert_search(
            input_path=config.search_results,
            search_engine=config.search_results_type,
            tmt_label=tmt_label,
            custom_mods=config.custom_to_unimod(),
            output_file=internal_search_file,
        )
        shutil.rmtree(config.output)
        pd.testing.assert_series_equal(search_results["MODIFIED_SEQUENCE"], exp_data["MODIFIED_SEQUENCE"])

    def test_no_mods_maxquant(self, exp_data_no_mods):
        """Test of no custom modifications in Maxquant results."""
        config = Config()
        config.read(Path(__file__).parent / "configs/no_custom_mods.json")

        msms_output = config.output / "msms"
        msms_output.mkdir(exist_ok=True, parents=True)
        internal_search_file = msms_output / "msms.prosit"
        tmt_label = config.tag
        search_results = pp.convert_search(
            input_path=config.search_results,
            search_engine=config.search_results_type,
            tmt_label=tmt_label,
            custom_mods=config.custom_to_unimod(),
            output_file=internal_search_file,
        )
        shutil.rmtree(config.output)
        pd.testing.assert_series_equal(search_results["MODIFIED_SEQUENCE"], exp_data_no_mods["MODIFIED_SEQUENCE"])


@pytest.fixture
def exp_data():
    """Simple data of modified sequences with new mods in internal format."""
    return pd.DataFrame(
        {
            "MODIFIED_SEQUENCE": [
                "AAAAAAYY",
                "AAAAPPWAC[UNIMOD:4]FAAV",
                "AAAAVVSGPKRGRKKP",
                "AAAAVVN[UNIMOD:7]SGPKRGRKKP",
                "AAAC[UNIMOD:4]RFVQ",
                "AAAC[UNIMOD:4]RFVQ",
                "AAAC[UNIMOD:4]RFVQ",
                "AAADQNNYHYL",
                "AAAGLVLIRFIFLVTLWVK",
                "AAAGRKVTLTTNLLLVG",
            ]
        }
    )


@pytest.fixture
def exp_data_no_mods():
    """Simple data of modified sequences in internal format."""
    return pd.DataFrame(
        {
            "MODIFIED_SEQUENCE": [
                "AAAAAAYY",
                "AAAAPPWAC[UNIMOD:4]FAAV",
                "AAAAVVSGPKRGRKKP",
                "AAAAVVNSGPKRGRKKP",
                "AAAC[UNIMOD:4]RFVQ",
                "AAAC[UNIMOD:4]RFVQ",
                "AAAC[UNIMOD:4]RFVQ",
                "AAADQNNYHYL",
                "AAAGLVLIRFIFLVTLWVK",
                "AAAGRKVTLTTNLLLVG",
            ]
        }
    )
