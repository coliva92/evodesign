from typing import Optional
from datetime import datetime
import os


class DirectoryManager:

    GIT_COMMIT_HASH = "a79ffafa989506f6cc3dfbd98be82d6942ec7908"

    def __init__(
        self,
        output_dir: str,
        jobname: Optional[str] = None,
    ):
        self._output_dir = os.path.abspath(output_dir)
        if jobname is None:
            today = datetime.today().strftime("%Y%m%d_%H%M%S")
            jobname = f"evodesign_{today}"
        self.jobname = jobname
        self.path = os.path.join(self._output_dir, self.jobname)
        self.results_npz_path = os.path.join(self.path, "results.npz")
        self.pymoo_algorithm_bin_path = os.path.join(self.path, "pymoo_algorithm.bin")
        self.initial_rng_state_path = os.path.join(self.path, "initial_rng_state.txt")
        self.last_rng_state_path = os.path.join(self.path, "last_rng_state.txt")
        self.prediction_pdbs_dir = os.path.join(self.path, "prediction_pdbs")
        self.settings_json_path = os.path.join(self.path, "settings.json")
        self.git_commit_hash_path = os.path.join(self.path, "git_commit_hash.txt")
        self.predictor_input_dir = os.path.join(self.path, "predictor_input")
        self.predictor_output_dir = os.path.join(self.path, "predictor_output")
        self.profile_path = os.path.join(self.path, "aa_profile.txt")
        self._essential_files = {
            self.results_npz_path,
            self.pymoo_algorithm_bin_path,
            self.initial_rng_state_path,
            self.last_rng_state_path,
            self.settings_json_path,
            self.git_commit_hash_path,
            self.profile_path
        }
        return

    def is_essential_file_or_folder(
        self,
        filepath: str,
    ) -> bool:
        return filepath in self._essential_files
