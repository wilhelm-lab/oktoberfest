import unittest
from pathlib import Path

from oktoberfest.utils import JobPool, ProcessStep


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
