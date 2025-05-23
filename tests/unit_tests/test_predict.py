import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, create_autospec, patch

import numpy as np
import pandas as pd
from numpy.testing import assert_almost_equal

import oktoberfest
from oktoberfest.data import Spectra
from oktoberfest.data.spectra import FragmentType
from oktoberfest.predict import Koina, Predictor
from oktoberfest.utils import Config

DATA_PATH = Path(__file__).parent / "data"


class TestTMTProsit(unittest.TestCase):
    """Test class for TMT model predictions."""

    def test_prosit_tmt(self):
        """Test retrieval of predictions from prosit tmt models via koina."""
        meta_df = pd.read_csv(Path(__file__).parent / "data" / "predictions" / "library_input.csv")
        var = Spectra._gen_vars_df()
        library = Spectra(obs=meta_df, var=var)
        library.strings_to_categoricals()

        intensity_predictor = Predictor.from_koina(
            model_name="Prosit_2020_intensity_TMT",
            server_url="koina.wilhelmlab.org:443",
            ssl=True,
        )
        intensity_predictor.predict_intensities(data=library)

        irt_predictor = Predictor.from_koina(
            model_name="Prosit_2020_irt_TMT",
            server_url="koina.wilhelmlab.org:443",
            ssl=True,
        )
        irt_predictor.predict_rt(data=library)

        library_expected = Spectra.from_hdf5(Path(__file__).parent / "data" / "predictions" / "library_output.h5ad.gz")

        assert_almost_equal(
            library.get_matrix(FragmentType.PRED).toarray(),
            library_expected.get_matrix(FragmentType.PRED).toarray(),
            decimal=6,
        )
        # explicitly set this to int64 as the default on windows is int32
        # which causes the type check to fail on windows
        library.var["num"] = library.var["num"].astype("int64")
        library.var["charge"] = library.var["charge"].astype("int64")

        pd.testing.assert_frame_equal(library.obs, library_expected.obs)
        pd.testing.assert_frame_equal(library.var, library_expected.var)


# TODO test kwarg passing
class TestPredictorBehavioral(unittest.TestCase):
    """Behavioral tests of Predictor class."""

    @classmethod
    def setUpClass(cls):  # noqa: 402
        cls.model_name = "Prosit_2019_intensity"
        cls.model_type = "intensity"

        cls.temp_dir = Path(tempfile.mkdtemp())
        cls.data_dir = cls.temp_dir / "data"
        cls.data_dir.mkdir()

        cls.mock_config = create_autospec(Config, instance=True)
        cls.mock_config.data = {}
        cls.mock_config.data["models"] = {cls.model_type: cls.model_name}
        cls.mock_config.output = cls.temp_dir

        cls.mock_koina = create_autospec(Koina, instance=False)
        cls.mock_spectra = create_autospec(Spectra, instance=True)
        cls.intensities = np.array([[0.0, 0.0, -1.0], [1.0, 0, -1.0], [1.0, 0.0, 0.0]])
        cls.ion_annotations = np.array(
            [["y1+1", "y1+2", "y1+3"], ["y1+1", "y1+2", "y1+3"], ["y1+1", "y1+2", "y1+3"]], dtype=object
        )
        cls.retention_times = np.array([30.0, 100.0, 160.0, 140.0, -2.0, 17.0])
        cls.chunk_idx = [pd.Index([0, 1, 2]), pd.Index([3, 4, 5])]
        cls.ce_range = (19, 50)

    @classmethod
    def tearDownClass(cls):  # noqa: 402
        shutil.rmtree(cls.temp_dir)

    @patch("oktoberfest.predict.predictor.Koina")
    def test_from_koina(self, mock_koina):
        """Test Koina constructor for Predictor."""
        Predictor.from_koina(model_name=self.model_name)
        assert mock_koina is oktoberfest.predict.predictor.Koina
        mock_koina.assert_called_once()

    @patch("oktoberfest.predict.predictor.DLomix")
    def test_from_dlomix(self, mock_dlomix):
        """Test DLomix constructor for Predictor."""
        if mock_dlomix is None:
            self.assertTrue(True)
            return
        Predictor.from_dlomix(
            model_type=self.model_type,
            model_path=self.temp_dir / "prosit_baseline.keras",
            output_path=self.temp_dir / "dlomix_output",
            batch_size=1024,
        )
        assert mock_dlomix is oktoberfest.predict.predictor.DLomix
        mock_dlomix.assert_called_once()

    @patch("oktoberfest.predict.predictor.Koina")
    def test_koina_from_config(self, mock_koina):
        """Test config constructor for Predictor with Koina."""
        self.mock_config.predict_intensity_locally = False
        Predictor.from_config(self.mock_config, model_type=self.model_type)
        assert mock_koina is oktoberfest.predict.predictor.Koina
        mock_koina.assert_called_once()

    @patch("oktoberfest.predict.predictor.DLomix")
    def test_dlomix_from_config(self, mock_dlomix):
        """Test config constructor for Predictor with DLomix."""
        if mock_dlomix is None:
            self.assertTrue(True)
            return
        self.mock_config.predict_intensity_locally = True
        self.mock_config.download_baseline_intensity_predictor = False
        self.mock_config.dlomix_inference_batch_size = 1024
        Predictor.from_config(self.mock_config, model_type=self.model_type)
        assert mock_dlomix is oktoberfest.predict.predictor.DLomix
        mock_dlomix.assert_called_once()

    @patch("oktoberfest.predict.predictor.DLomix")
    def test_download_new_model(self, mock_dlomix):
        """Test if new baseline model is downloaded if requested."""
        if mock_dlomix is None:
            self.assertTrue(True)
            return
        self.mock_config.download_baseline_intensity_predictor = True
        Predictor.from_config(self.mock_config, model_type=self.model_type)
        assert mock_dlomix is oktoberfest.predict.predictor.DLomix
        mock_dlomix.assert_called_once_with(
            model_type=self.model_type,
            model_path=self.data_dir / "dlomix/prosit_baseline_model.keras",
            output_path=self.data_dir / "dlomix",
            batch_size=self.mock_config.dlomix_inference_batch_size,
            download=True,
        )

    def test_predict_intensities_at_once(self):
        """Test if predict_intensities does the right steps when chunk_idx=None."""
        # TODO add state-based test
        predictor = Predictor(self.mock_koina, model_name=self.model_name)
        predictor._predictor.predict = MagicMock(
            return_value={"intensities": self.intensities, "annotation": self.ion_annotations}
        )
        predictor.predict_intensities(self.mock_spectra)
        predictor._predictor.predict.assert_called_once_with(data=self.mock_spectra)
        self.mock_spectra.add_intensities.assert_called_once_with(
            self.intensities, self.ion_annotations, fragment_type=FragmentType.PRED
        )

    def test_predict_intensities_chunked(self):
        """Test if predict_intensities does the right steps with chunk_idx."""
        # TODO add state-based test
        predictor = Predictor(self.mock_koina, model_name=self.model_name)
        predictor.predict_in_chunks = MagicMock(
            return_value={
                "intensities": [self.intensities, self.intensities],
                "annotation": [self.ion_annotations, self.ion_annotations],
            }
        )
        predictor.predict_intensities(self.mock_spectra, chunk_idx=self.chunk_idx)
        predictor.predict_in_chunks.assert_called_once_with(data=self.mock_spectra, chunk_idx=self.chunk_idx, xl=False)
        self.mock_spectra.add_list_of_predicted_intensities.assert_called_once_with(
            [self.intensities, self.intensities], [self.ion_annotations, self.ion_annotations], self.chunk_idx
        )

    def test_predict_rt(self):
        """Test iRT prediction."""
        # TODO add state-based test
        predictor = Predictor(self.mock_koina, model_name=self.model_name)
        predictor._predictor.predict = MagicMock(return_value={"irt": self.retention_times})
        predictor.predict_rt(self.mock_spectra)
        predictor._predictor.predict.assert_called_once_with(self.mock_spectra)
        self.mock_spectra.add_column.assert_called_once_with(self.retention_times, name="PREDICTED_IRT")

    def test_predict_at_once(self):
        """Test prediction in one go."""
        predictor = Predictor(self.mock_koina, model_name=self.model_name)
        result = {"intensities": self.intensities, "annotation": self.ion_annotations}
        predictor._predictor.predict = MagicMock(return_value=result)
        output = predictor.predict_at_once(self.mock_spectra)
        predictor._predictor.predict.assert_called_once_with(self.mock_spectra)
        self.assertEqual(output, result)

    def test_predict_in_chunks(self):
        """Test prediction in chunks."""
        # TODO add state-based test
        predictor = Predictor(self.mock_koina, model_name=self.model_name)
        predictor._predictor.predict = MagicMock(
            side_effect=[
                {"intensities": self.intensities, "annotation": self.ion_annotations},
                {"intensities": self.intensities, "annotation": self.ion_annotations},
            ]
        )
        output = predictor.predict_in_chunks(self.mock_spectra, chunk_idx=self.chunk_idx)
        self.assertEqual(
            output,
            {
                "intensities": [self.intensities, self.intensities],
                "annotation": [self.ion_annotations, self.ion_annotations],
            },
        )

    @patch("oktoberfest.predict.predictor._prepare_alignment_df")
    @patch("oktoberfest.predict.predictor._alignment")
    def test_ce_calibration(self, mock_alignment, mock_prepare_alignment_df):
        """Test CE calibration."""
        # TODO add state-based test
        mock_prepare_alignment_df.return_value = self.mock_spectra
        self.mock_spectra.obs = pd.DataFrame({"PEPTIDE_LENGTH": [9, 11, 8, 9, 11, 8]})
        predictor = Predictor(self.mock_koina, model_name=self.model_name)
        predictor.predict_intensities = MagicMock()
        alignment_library = predictor.ce_calibration(
            library=self.mock_spectra, ce_range=self.ce_range, group_by_charge=False, model_name="prosit"
        )
        mock_prepare_alignment_df.assert_called_once_with(
            self.mock_spectra, ce_range=self.ce_range, group_by_charge=False, xl=False
        )
        predictor.predict_intensities.assert_called_once_with(
            data=self.mock_spectra, chunk_idx=None, keep_dataset=False, model_name="prosit", xl=False
        )
        mock_alignment.assert_called_once_with(self.mock_spectra, xl=False)
        self.assertEqual(alignment_library, self.mock_spectra)

    @patch("oktoberfest.predict.predictor.group_iterator", return_value="chunk_idx_dummy")
    @patch("oktoberfest.predict.predictor._prepare_alignment_df")
    @patch("oktoberfest.predict.predictor._alignment")
    def test_ce_calibration_alphapept(self, mock_alignment, mock_prepare_alignment_df, mock_group_iterator):
        """Test CE calibration with alphapept predictor."""
        mock_prepare_alignment_df.return_value = self.mock_spectra
        self.mock_spectra.obs = pd.DataFrame({"PEPTIDE_LENGTH": [9, 11, 8, 9, 11, 8]})
        predictor = Predictor(self.mock_koina, model_name="AlphaPept_ms2_generic")
        predictor.predict_intensities = MagicMock()
        alignment_library = predictor.ce_calibration(
            library=self.mock_spectra, ce_range=self.ce_range, group_by_charge=False
        )
        mock_prepare_alignment_df.assert_called_once_with(
            self.mock_spectra, ce_range=self.ce_range, group_by_charge=False, xl=False
        )
        predictor.predict_intensities.assert_called_once_with(
            data=self.mock_spectra, chunk_idx=list("chunk_idx_dummy"), keep_dataset=False, xl=False
        )
        mock_alignment.assert_called_once_with(self.mock_spectra, xl=False)
        self.assertEqual(alignment_library, self.mock_spectra)


class TestPredictorStateBased(unittest.TestCase):
    """State-based tests of Predictor class."""

    @classmethod
    def setUpClass(cls):  # noqa: D102
        cls.ce_range = (19, 50)
        cls.spectra = Spectra.from_hdf5(DATA_PATH / "spectra/test_spectra.hdf5")
        cls.irt_predictor = Predictor.from_koina(model_name="Prosit_2019_irt")
        cls.predictor = Predictor.from_koina(model_name="Prosit_2019_intensity")

    def test_ce_calibration(self):
        """State-based CE calibration test."""
        # TODO
        pass
        """
        alignment_library = self.predictor.ce_calibration(
            library=self.spectra, ce_range=self.ce_range, group_by_charge=False
        )
        """

    def test_ce_calibration_chunked(self):
        """State-based CE calibration test."""
        # TODO
        pass
        """
        alignment_library = self.predictor.ce_calibration(
            library=self.spectra, ce_range=self.ce_range, group_by_charge=True
        )
        """


class TestLocalPrediction(unittest.TestCase):
    """Test class for local prediction."""

    @classmethod
    def setUpClass(cls):  # noqa: 402
        cls.model_name = "Prosit_2019_intensity"
        cls.model_type = "intensity"

        cls.temp_dir = Path(tempfile.mkdtemp())
        cls.data_dir = cls.temp_dir / "data"
        cls.data_dir.mkdir()

        cls.mock_config = create_autospec(Config, instance=True)
        cls.mock_config.data = {}
        cls.mock_config.data["models"] = {cls.model_type: cls.model_name}
        cls.mock_config.output = cls.temp_dir

        cls.mock_spectra = create_autospec(Spectra, instance=True)
        cls.intensities = np.array([[0.0, 0.0, -1.0], [1.0, 0, -1.0], [1.0, 0.0, 0.0]])
        cls.ion_annotations = np.array(
            [["y1+1", "y1+2", "y1+3"], ["y1+1", "y1+2", "y1+3"], ["y1+1", "y1+2", "y1+3"]], dtype=object
        )
        cls.retention_times = np.array([30.0, 100.0, 160.0, 140.0, -2.0, 17.0])
        cls.chunk_idx = [pd.Index([0, 1, 2]), pd.Index([3, 4, 5])]
        cls.ce_range = (25, 30)

    @patch("oktoberfest.predict.predictor.DLomix")
    def test_predict_rt(self, mock_dlomix):
        """Test iRT prediction."""
        # TODO add state-based test
        if mock_dlomix is None:
            self.assertTrue(True)
            return
        predictor = Predictor(mock_dlomix, model_name=self.model_name)
        predictor._predictor.predict = MagicMock(return_value={"irt": self.retention_times})
        predictor.predict_rt(self.mock_spectra)
        predictor._predictor.predict.assert_called_once_with(self.mock_spectra)
        self.mock_spectra.add_column.assert_called_once_with(self.retention_times, name="PREDICTED_IRT")


class TestRefinementLearning(unittest.TestCase):
    """Test class for refinement learning."""
