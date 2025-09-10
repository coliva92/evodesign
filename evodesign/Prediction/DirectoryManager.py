import os
import shutil


class DirectoryManager:

    def __init__(
        self,
        prediction_pdbs_dir: str,
        model_input_dir: str,
        model_output_dir: str,
        prefix: str = "tmp_prediction",
    ):
        self.prediction_pdbs_dir = prediction_pdbs_dir
        self.model_input_dir = model_input_dir
        self.model_output_dir = model_output_dir
        self.prefix = prefix

    def _empty_folder(self, folder_path: str) -> None:
        if not os.path.isdir(folder_path):
            return
        for filename in os.listdir(folder_path):
            file_path = os.path.join(folder_path, filename)
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
                continue
            # if not file or symlink, then it's a folder
            shutil.rmtree(file_path)

    def empty_folders_content(self) -> None:
        self._empty_folder(self.prediction_pdbs_dir)
        self._empty_folder(self.model_input_dir)
        self._empty_folder(self.model_output_dir)
