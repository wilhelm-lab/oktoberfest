import tempfile
import unittest
import warnings
from pathlib import Path

import anndata
import numpy as np
import pandas as pd

from oktoberfest import preprocessing as pp
from oktoberfest.data.spectra import FragmentType, Spectra
from oktoberfest.plotting.plotting import joint_plot


class TestLibraryObs(unittest.TestCase):
    """Test the obs of spectral libraries get a string index and are no views (no AnnData warnings)."""

    def test_gen_lib(self):
        """Test the library read from a peptide list has a string index."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            peptides = Path(tmp_dir) / "peptides.csv"
            pd.DataFrame(
                {
                    "modified_sequence": ["KIIDRAITSL", "KIEKLKVEL", "ASPTQPIQL"],
                    "collision_energy": [34, 35, 34],
                    "precursor_charge": [2, 2, 1],
                }
            ).to_csv(peptides, sep="\t", index=False)
            with warnings.catch_warnings():
                warnings.simplefilter("error", anndata.ImplicitModificationWarning)
                library = pp.gen_lib(peptides)

        self.assertListEqual(list(library.obs_names), ["0", "1", "2"])
        self.assertListEqual(library.obs["PROTEINS"].tolist(), ["unknown"] * 3)
        self.assertEqual(library.n_vars, len(Spectra._gen_vars_df()))

    def test_process_and_filter_spectra_data(self):
        """Test peptides too long for the model are removed and the result is a copy with masses."""
        obs = pd.DataFrame(
            {
                "MODIFIED_SEQUENCE": ["AAACRFVQ", "RM[UNIMOD:35]PCHKPYL", "A" * 40],
                "PRECURSOR_CHARGE": [1, 2, 2],
            },
            index=["0", "1", "2"],
        )
        library = Spectra(obs=obs, var=Spectra._gen_vars_df())
        with warnings.catch_warnings():
            warnings.simplefilter("error", anndata.ImplicitModificationWarning)
            library = pp.process_and_filter_spectra_data(library, "Prosit_2020_intensity_HCD")

        self.assertFalse(library.is_view)
        self.assertListEqual(
            library.obs["MODIFIED_SEQUENCE"].tolist(), ["AAAC[UNIMOD:4]RFVQ", "RM[UNIMOD:35]PC[UNIMOD:4]HKPYL"]
        )
        self.assertListEqual(library.obs["SEQUENCE"].tolist(), ["AAACRFVQ", "RMPCHKPYL"])
        self.assertTrue((library.obs["MASS"] > 0).all())

    def test_annotate_spectral_library(self):
        """Test annotating PSMs with a default integer index gives a string obs index and the matched peaks."""
        psms = pd.DataFrame(
            {
                "RAW_FILE": ["File1"],
                "SCAN_NUMBER": [1],
                "MODIFIED_SEQUENCE": ["PEPTIDEK"],
                "SEQUENCE": ["PEPTIDEK"],
                "PRECURSOR_CHARGE": [2],
                "PEPTIDE_LENGTH": [8],
                "MASS_ANALYZER": ["FTMS"],
                "FRAGMENTATION": ["HCD"],
                "INTENSITIES": [np.array([0.2, 1.0, 0.5])],
                "MZ": [np.array([147.11280, 227.10263, 276.15539])],
            }
        )
        with warnings.catch_warnings():
            warnings.simplefilter("error", anndata.ImplicitModificationWarning)
            library = pp.annotate_spectral_library(psms, mass_tol=20, unit_mass_tol="ppm")

        self.assertListEqual(list(library.obs_names), ["0"])
        self.assertNotIn("INTENSITIES", library.obs.columns)
        # possible but unobserved fragments get a placeholder intensity of 1e-7
        intensities = library.get_matrix(FragmentType.RAW).toarray()[0]
        self.assertListEqual(library.var_names[intensities > 1e-3].tolist(), ["y1+1", "y2+1", "b2+1"])


class TestRemoveDuplicates(unittest.TestCase):
    """Test removing (peptide, charge, CE) duplicates."""

    def test_keeps_first_n(self):
        """Test the first n spectra of each (peptide, charge, CE) are kept and the helper column is dropped."""
        obs = pd.DataFrame(
            {
                "SEQUENCE": ["PEPTIDEK"] * 4 + ["LESLIEK"] * 2,
                "PRECURSOR_CHARGE": [2, 2, 2, 3, 2, 2],
                "COLLISION_ENERGY": [30, 30, 30, 30, 30, 35],
                "SCAN_NUMBER": [1, 2, 3, 4, 5, 6],
            },
            index=[str(i) for i in range(6)],
        )
        spectra = Spectra(obs=obs, var=Spectra._gen_vars_df())
        spectra.remove_duplicates(2)

        self.assertListEqual(spectra.obs["SCAN_NUMBER"].tolist(), [1, 2, 4, 5, 6])
        self.assertNotIn("duplicate_count", spectra.obs.columns)


class TestJointPlot(unittest.TestCase):
    """Test the joint plot of rescoring with and without peptide property predictions."""

    def test_partial_overlap(self):
        """Test PSMs found by only one of the two rescorings (NaN scores after the outer merge) are dropped."""
        target = pd.DataFrame(
            {
                "PSMId": ["F1-15-TAIASPEK-1-5", "F2-59-LGLTKLQLH-3-9", "F1-5-YLNPLRTK-1-5"],
                "q-value": [0.001, 0.002, 0.003],
                "score": [0.4, 0.5, 0.6],
                "peptide": ["TAIASPEK", "LGLTKLQLH", "YLNPLRTK"],
            }
        )
        decoy = pd.DataFrame(
            {
                "PSMId": ["F1-11-KLYNANYIK-3-7", "F1-24-EFAVEVLK-2-4"],
                "q-value": [0.006, 0.003],
                "score": [-0.1, -0.5],
                "peptide": ["KLYNANYIK", "EFAVEVLK"],
            }
        )
        andromeda_target = target.iloc[:2].assign(score=[0.2, 0.3])
        andromeda_decoy = decoy.iloc[1:].assign(score=[-0.4])
        with tempfile.TemporaryDirectory() as tmp_dir:
            filename = Path(tmp_dir) / "joint_plot.png"
            joint_plot(target, decoy, andromeda_target, andromeda_decoy, level="psm", filename=filename)
            self.assertTrue(filename.exists())
