import logging
import os

logger = logging.getLogger(__name__)


class ProcessStep:
    """Init a ProcessStep object given raw file path and step name."""

    def __init__(self, out_path: str, step_name: str):
        """
        Init raw file path and step name.

        :param raw_path: path to raw file
        :param step_name: name of the current step
        """
        self.out_path = out_path
        self.step_name = step_name

    def _get_proc_folder_path(self) -> str:
        """Get proc folder path."""
        return os.path.join(self.out_path, "proc")

    def _get_done_file_path(self) -> str:
        """Get done if a file is done."""
        return os.path.join(self._get_proc_folder_path(), self.step_name + ".done")

    def is_done(self) -> bool:
        """Retrun True if file is done."""
        if not os.path.isdir(self._get_proc_folder_path()):
            os.makedirs(self._get_proc_folder_path())

        if os.path.isfile(self._get_done_file_path()):
            logger.info(f"Skipping {self.step_name} step because {self._get_done_file_path()} was found.")
            return True
        else:
            return False

    def mark_done(self):
        """Mark file as done."""
        open(self._get_done_file_path(), "w").close()
