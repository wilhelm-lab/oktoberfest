import json
import shutil
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from oktoberfest.utils import Config, JobPool, ProcessStep
from oktoberfest.utils.quantification import apply_quant


def add_one(i: int):
    """Test function for multiprocessing pool."""
    return i + 1


class TestJobPool(unittest.TestCase):
    """Test the JobPool class."""

    def test_jobpool(self):
        """Unit test for starting and joining multiprocessing pool."""
        pool = JobPool(2)
        for i in range(5):
            pool.apply_async(add_one, [i])
        self.assertEqual(pool.check_pool(), [1, 2, 3, 4, 5])


class TestProcessStep(unittest.TestCase):
    """Test the JobPool class."""

    def test_process_step(self):
        """Unit test for starting and joining multiprocessing pool."""
        proc_step = ProcessStep(out_path=str(Path(__file__).parent), step_name="test_step")
        proc_step = ProcessStep(out_path=Path(__file__).parent, step_name="test_step")

        self.assertFalse(proc_step.is_done())
        proc_step.mark_done()
        self.assertTrue(proc_step.is_done())
        proc_step_file = proc_step._get_done_file_path()
        self.assertTrue(proc_step_file.is_file())
        proc_step_file.unlink()


class TestConfig(unittest.TestCase):
    """Test the Config class."""

    @classmethod
    def setUpClass(cls):  # noqa: D102
        cls.config_path = Path(__file__).parents[1] / "configs/rescoring_local_prediction.json"
        cls.temp_dir = Path(tempfile.mkdtemp())

    @classmethod
    def tearDownClass(cls):  # noqa: D102
        shutil.rmtree(cls.temp_dir)

    def test_check_model_availability(self):
        """Test if invalid model is being checked."""
        with open(self.config_path) as f:
            raw_config = json.load(f)
        raw_config["models"]["intensity"] = "garbage"
        garbage_config_file = self.temp_dir / "garbage_config.json"
        with open(garbage_config_file, "w+") as f:
            json.dump(raw_config, f)

        conf = Config()
        conf.read(garbage_config_file)
        with self.assertRaises(ValueError):
            conf.check()
        garbage_config_file.unlink()

    def test_num_threads_defaults_to_one(self):
        """Test that num_threads defaults to one when unset."""
        conf = Config()
        conf.data = {}
        self.assertEqual(conf.num_threads, 1)

    def test_num_threads_uses_config_key(self):
        """Test that num_threads uses the numThreads config key."""
        conf = Config()
        conf.data = {"numThreads": 4}
        self.assertEqual(conf.num_threads, 4)


class TestQuant(unittest.TestCase):
    """Test the quantification done by calling picked-group-fdr."""

    def test_picked_group_fdr_maxquant(self):
        """Testing picked_group_fdr quantification with msfragger search results."""
        config = Config()
        config.data = {
            "inputs": {
                "search_results": Path("./data/quantification/mq"),
                "search_results_type": "maxquant",
                "library_input": Path("./data/quantification/example.fasta"),
            },
            "output": Path("./data/quantification"),
            "fdr_estimation_method": "percolator",
            "fastaDigestOptions": {
                "digestion": "full",
                "missedCleavages": 2,
                "minLength": 7,
                "maxLength": 60,
                "enzyme": "asp-n",
                "specialAas": "D",
                "db": "target",
            },
        }
        config.base_path = Path(__file__).parents[1]
        apply_quant(config)
        compare = pd.read_csv(Path(__file__).parents[1] / "data/quantification/mq_proteinGroups.txt", sep="\t")
        results = pd.read_csv(
            Path(__file__).parents[1] / "data/quantification/picked_group_fdr/rescore.proteinGroups.txt", sep="\t"
        )
        pd.testing.assert_frame_equal(results, compare)

    def test_picked_group_fdr_sage(self):
        """Testing picked_group_fdr quantification with sage search results."""
        config = Config()
        config.data = {
            "inputs": {
                "search_results": Path("./data/quantification/sage"),
                "search_results_type": "sage",
                "library_input": Path("./data/quantification/iprg2016_with_labels.fasta"),
            },
            "output": Path("./data/quantification"),
            "fdr_estimation_method": "percolator",
            "fastaDigestOptions": {
                "digestion": "full",
                "missedCleavages": 2,
                "minLength": 7,
                "maxLength": 60,
                "enzyme": "asp-n",
                "specialAas": "D",
                "db": "target",
            },
        }
        config.base_path = Path(__file__).parents[1]
        # TODO add data for testing
        # apply_quant(config)

    def test_picked_group_fdr_fragpipe(self):
        """Testing picked_group_fdr quantification with msfragger search results."""
        config = Config()
        config.data = {
            "inputs": {
                "search_results": Path("./data/quantification/fragpipe"),
                "search_results_type": "msfragger",
                "library_input": Path("./data/quantification/iprg2016_with_labels.fasta"),
            },
            "output": Path("./data/quantification"),
            "fdr_estimation_method": "percolator",
            "fastaDigestOptions": {
                "digestion": "full",
                "missedCleavages": 2,
                "minLength": 7,
                "maxLength": 60,
                "enzyme": "asp-n",
                "specialAas": "D",
                "db": "target",
            },
        }
        config.base_path = Path(__file__).parents[1]
        # TODO add data for testing
        # apply_quant(config)
