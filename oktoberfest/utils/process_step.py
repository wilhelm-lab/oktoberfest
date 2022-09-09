import logging
import os

logger = logging.getLogger(__name__)


class ProcessStep:
    def __init__(self, out_path, step_name):
        self.out_path = out_path
        self.step_name = step_name

    def _get_proc_folder_path(self):
        return os.path.join(self.out_path, "proc")

    def _get_done_file_path(self):
        return os.path.join(self._get_proc_folder_path(), self.step_name + ".done")

    def is_done(self):
        if not os.path.isdir(self._get_proc_folder_path()):
            os.makedirs(self._get_proc_folder_path())

        if os.path.isfile(self._get_done_file_path()):
            logger.info(f"Skipping {self.step_name} step because {self._get_done_file_path()} was found.")
            return True
        else:
            return False

    def mark_done(self):
        open(self._get_done_file_path(), "w").close()
