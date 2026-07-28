import os
from dataclasses import dataclass
from datetime import datetime
from typing import ClassVar, FrozenSet, Optional


@dataclass
class PathsContainer:

    PROJECT_VERSION: ClassVar[str] = "20260727184718"
    jobname: str
    path: str
    results_npz_path: str
    pymoo_algorithm_bin_path: str
    initial_rng_state_path: str
    last_rng_state_path: str
    prediction_pdbs_dir: str
    settings_json_path: str
    project_version_path: str
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
        version = os.path.join(base_path, "version.txt")
        profile = os.path.join(base_path, "aa_profile.txt")
        essential_files = frozenset(
            {results, pymoo_bin, init_rng, last_rng, settings, version, profile}
        )
        return cls(
            jobname=jobname,
            path=base_path,
            results_npz_path=results,
            pymoo_algorithm_bin_path=pymoo_bin,
            initial_rng_state_path=init_rng,
            last_rng_state_path=last_rng,
            prediction_pdbs_dir=os.path.join(base_path, "prediction_pdbs"),
            settings_json_path=settings,
            project_version_path=version,
            predictor_input_dir=os.path.join(base_path, "predictor_input"),
            predictor_output_dir=os.path.join(base_path, "predictor_output"),
            profile_path=profile,
            essential_files=essential_files,
        )
