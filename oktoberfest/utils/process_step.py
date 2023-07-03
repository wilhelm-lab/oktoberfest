import logging
from pathlib import Path
from typing import Union

logger = logging.getLogger(__name__)


class ProcessStep:
    """Init a ProcessStep object given raw file path and step name."""

    def __init__(self, out_path: Union[str, Path], step_name: str):
        """
        Init raw file path and step name.

        :param out_path: path to raw file
        :param step_name: name of the current step
        """
        if isinstance(out_path, str):
            out_path = Path(out_path)
        self.out_path = out_path
        self.step_name = step_name

    def _get_proc_folder_path(self) -> Path:
        """Get proc folder path."""
        return self.out_path / "proc"

    def _get_done_file_path(self) -> Path:
        """Get done if a file is done."""
        return self._get_proc_folder_path() / f"{self.step_name}.done"

    def is_done(self) -> bool:
        """Retrun True if file is done."""
        self._get_proc_folder_path().mkdir(exist_ok=True)

        if self._get_done_file_path().is_file():
            logger.info(f"Skipping {self.step_name} step because {self._get_done_file_path()} was found.")
            return True
        else:
            return False

    def mark_done(self):
        """Mark file as done."""
        open(self._get_done_file_path(), "w").close()
