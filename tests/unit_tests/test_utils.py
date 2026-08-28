import ast
import inspect
import json
import shutil
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from oktoberfest import preprocessing as pp
from oktoberfest import runner
from oktoberfest import rescore as re
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
        pool.check_pool()


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


class TestFeatureAndMatcherOptions(unittest.TestCase):
    """The config options this branch adds, and the fact that they actually reach the pipeline."""

    @staticmethod
    def _config(data):
        config = Config()
        config.data = data
        return config

    def test_sc_features_off_by_default(self):
        """The single-cell feature families must not change the default feature set."""
        self.assertFalse(self._config({}).sc_features)
        self.assertFalse(self._config({"sc_features": False}).sc_features)

    def test_sc_features_can_be_switched_on(self):
        """The explicit flag switches them on."""
        self.assertTrue(self._config({"sc_features": True}).sc_features)

    def test_all_features_implies_sc_features(self):
        """The allFeatures flag means all features, so it must not skip the sc_feature families."""
        self.assertTrue(self._config({"allFeatures": True}).sc_features)
        self.assertTrue(self._config({"allFeatures": True, "sc_features": False}).sc_features)

    @staticmethod
    def _forwarded_kwargs(func, callee: str) -> set:
        """Keyword arguments ``func`` passes to every call of ``callee`` in its own body.

        This is a wiring assertion, not a behaviour one: the runner's feature-generation path needs a
        full prediction run to execute, so the regression this guards against -- an option that is read
        from the config but never forwarded, and therefore silently a no-op -- is checked statically.
        """
        tree = ast.parse(inspect.getsource(func))
        names = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                target = node.func
                name = target.attr if isinstance(target, ast.Attribute) else getattr(target, "id", None)
                if name == callee:
                    names |= {kw.arg for kw in node.keywords if kw.arg}
        return names

    def test_sc_features_reach_generate_features(self):
        """An option nobody forwards is a no-op; check it lands on rescore.generate_features."""
        config = self._config({"sc_features": True})
        self.assertTrue(config.sc_features)
        signature = inspect.signature(re.generate_features)
        self.assertIn("sc_features", signature.parameters)
        self.assertIs(signature.parameters["sc_features"].default, False)

    def test_the_runner_forwards_sc_features(self):
        """Both the original and the rescore feature sets must get the flag."""
        source = inspect.getsource(runner._calculate_features)
        self.assertEqual(source.count("sc_features=config.sc_features"), 2)
        self.assertIn("sc_features", self._forwarded_kwargs(runner._calculate_features, "generate_features"))

    def test_the_runner_forwards_the_matching_options(self):
        """matching_method_params was once read from the config but never forwarded; guard both."""
        forwarded = self._forwarded_kwargs(runner._annotate_and_get_library, "annotate_spectral_library")
        self.assertIn("matching_method", forwarded)
        self.assertIn("matching_method_params", forwarded)

    def test_matching_method_defaults_to_nearest(self):
        """The default resolver reproduces the legacy closest-m/z behaviour."""
        self.assertEqual(self._config({}).matching_method, "nearest")
        self.assertEqual(self._config({}).matching_method_params, {})

    def test_matching_method_is_read_from_the_config(self):
        """Both the resolver name and its keyword arguments are read through."""
        config = self._config({"matching_method": "global_ransac", "matching_method_params": {"unique_peak": False}})
        self.assertEqual(config.matching_method, "global_ransac")
        self.assertEqual(config.matching_method_params, {"unique_peak": False})

    def test_unknown_matching_method_is_rejected_before_the_run(self):
        """A typo must fail at config check, not hours later inside annotation."""
        config = self._config({"matching_method": "nearest_neighbour"})
        with self.assertRaises(ValueError) as ctx:
            config.check()
        self.assertIn("nearest_neighbour", str(ctx.exception))

    def test_matching_method_reaches_annotation(self):
        """Both options must be forwarded by the annotation entry point, or they are a no-op."""
        signature = inspect.signature(pp.annotate_spectral_library)
        self.assertIn("matching_method", signature.parameters)
        self.assertIn("matching_method_params", signature.parameters)
        self.assertEqual(signature.parameters["matching_method"].default, "nearest")


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
