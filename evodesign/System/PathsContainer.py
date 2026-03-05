from dataclasses import dataclass
from typing import Optional, ClassVar, FrozenSet
from datetime import datetime
import os


@dataclass
class PathsContainer:

    GIT_COMMIT_HASH: ClassVar[str] = "a79ffafa989506f6cc3dfbd98be82d6942ec7908"
    jobname: str
    path: str
    results_npz_path: str
    pymoo_algorithm_bin_path: str
    initial_rng_state_path: str
    last_rng_state_path: str
    prediction_pdbs_dir: str
    settings_json_path: str
    git_commit_hash_path: str
    predictor_input_dir: str
    predictor_output_dir: str
    profile_path: str
    essential_files: FrozenSet[str]

    @classmethod
    def create(cls, output_dir: str, jobname: Optional[str] = None) -> "PathsContainer":
        abs_output_dir = os.path.abspath(output_dir)
        if jobname is None:
            today = datetime.today().strftime("%Y%m%d_%H%M%S")
            jobname = f"evodesign_{today}"
        base_path = os.path.join(abs_output_dir, jobname)
        results = os.path.join(base_path, "results.npz")
        pymoo_bin = os.path.join(base_path, "pymoo_algorithm.bin")
        init_rng = os.path.join(base_path, "initial_rng_state.txt")
        last_rng = os.path.join(base_path, "last_rng_state.txt")
        settings = os.path.join(base_path, "settings.json")
        git_hash = os.path.join(base_path, "git_commit_hash.txt")
        profile = os.path.join(base_path, "aa_profile.txt")
        essential_files = frozenset({
            results, pymoo_bin, init_rng, last_rng, settings, git_hash, profile
        })
        return cls(
            jobname=jobname,
            path=base_path,
            results_npz_path=results,
            pymoo_algorithm_bin_path=pymoo_bin,
            initial_rng_state_path=init_rng,
            last_rng_state_path=last_rng,
            prediction_pdbs_dir=os.path.join(base_path, "prediction_pdbs"),
            settings_json_path=settings,
            git_commit_hash_path=git_hash,
            predictor_input_dir=os.path.join(base_path, "predictor_input"),
            predictor_output_dir=os.path.join(base_path, "predictor_output"),
            profile_path=profile,
            essential_files=essential_files
        )
