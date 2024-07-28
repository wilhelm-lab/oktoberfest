from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata.tests.helpers import assert_equal

from oktoberfest.data.spectra import Spectra


class TestSpectra:
    """Test class for Spectra."""

    def test_read_and_write_hdf5(self, mini_spectra, tmpdir):
        """Test writing hdf5 files and reading them method."""
        spectra_path = Path(tmpdir / "spectra.hdf5")
        spec1 = mini_spectra
        spec1.write_as_hdf5(spectra_path)
        spec2 = Spectra.from_hdf5(spectra_path)
        assert_equal(spec1, spec2)

    def test_gen_vardf(self, var_df_cid, var_df_uvpd):
        """Test gen_vardf method."""
        var_df1 = Spectra._gen_vars_df()
        pd.testing.assert_frame_equal(var_df1, var_df_cid)
        var_df2 = Spectra._gen_vars_df(["x", "y", "z", "a", "b", "c"])
        pd.testing.assert_frame_equal(var_df2, var_df_uvpd)

    def test_assemble_for_parquet(self, mini_spectra, df_for_parquet):
        """Test assemlbe_for_parquet method."""
        df = mini_spectra.assemble_df_for_parquet()
        df = df.astype({"modified_sequence": "object"})
        df = df.astype({"modified_sequence": "category"})
        pd.testing.assert_frame_equal(
            df.reset_index(drop=True), df_for_parquet.reset_index(drop=True), check_column_type=False
        )


@pytest.fixture
def var_df_cid():
    """Returns var_df with y and b ions."""
    var_df = pd.read_csv(Path(__file__).parent / "data/spectra/vardf_CID.tsv", index_col="ion", sep="\t")
    return var_df


@pytest.fixture
def var_df_uvpd():
    """Returns var_df with z, y, x, c, b and a ions."""
    var_df = pd.read_csv(Path(__file__).parent / "data/spectra/vardf_UVPD.tsv", index_col="ion", sep="\t")
    return var_df


@pytest.fixture
def mini_spectra():
    """Returns spectra object with 10 entries."""
    s = Spectra.from_hdf5(Path(__file__).parent / "data/spectra/test_spectra.hdf5")
    return s


@pytest.fixture
def df_for_parquet():
    """Returns df with 10 enrtries in correct form for parquet/dlomix."""
    df = pd.read_csv(Path(__file__).parent / "data/spectra/df_for_parquet.csv", sep="\t", index_col="Unnamed: 0")
    df = df.astype({"method_nbr": "category", "modified_sequence": "category"})
    df["intensities_raw"] = df["intensities_raw"].map(
        lambda intens: np.fromstring(
            intens.replace("\n", "").replace("[", "").replace("]", "").replace("  ", " "), sep=" "
        )
    )
    df["precursor_charge_onehot"] = df["precursor_charge_onehot"].map(
        lambda intens: np.fromstring(
            intens.replace("\n", "").replace("[", "").replace("]", "").replace("  ", " "), dtype="int", sep=" "
        )
    )
    return df
