from typing import Optional
from datetime import datetime
import os


class WorkingFolder:

    GIT_COMMIT_HASH = "0a4e8c8d91311bba1b228d416eacb0c24d558a90"

    def __init__(self, output_dir: str, jobname: Optional[str] = None):
        self._output_dir = output_dir
        if jobname is None:
            today = datetime.today().strftime("%Y%m%d")
            jobname = f"evodesign_{today}"
        self.jobname = jobname
        self.path = os.path.join(self._output_dir, self.jobname)
        self.results_npz_path = os.path.join(self.path, "results.npz")
        self.pymoo_algorithm_bin_path = os.path.join(self.path, "pymoo_algorithm.bin")
        self.initial_rng_state_path = os.path.join(self.path, "initial_rng_state.txt")
        self.last_rng_state_path = os.path.join(self.path, "last_rng_state.txt")
        self.prediction_pdb_path = os.path.join(self.path, "prediction.pdb.tmp")
        self.settings_json_path = os.path.join(self.path, "settings.json")
        self.git_commit_hash_path = os.path.join(self.path, "git_commit_hash.txt")
