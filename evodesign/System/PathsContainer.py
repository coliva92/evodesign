import os
from dataclasses import dataclass
from datetime import datetime
from typing import ClassVar, FrozenSet, Optional




@dataclass
class PathsContainer:
    """
    Central registry for all file and directory paths used during a job's execution.
    This class ensures consistent path resolution across the system.
    
    Attributes:
        jobname (str): Unique identifier for the current run/job.
        path (str): The root directory for the job's outputs.
        results_npz_path (str): Path to the final numpy archive containing the results.
        pymoo_algorithm_bin_path (str): Path to the serialized pymoo optimizer state.
        initial_rng_state_path (str): Path to the starting random number generator state.
        latest_rng_state (str): Path to the latest random number generator state.
        prediction_pdbs_dir (str): Directory containing generated PDB structures.
        settings_json_path (str): Path to the serialized configuration for the run.
        program_version_path (str): Path to the file recording the program version.
        predictor_input_dir (str): Directory for staging inputs to the predictor model.
        predictor_output_dir (str): Directory for capturing outputs from the predictor 
            model.
        profile_path (str): Path to the probabilistic amino acid profile definition.
        essential_files (FrozenSet[str]): Set of critical output files that must remain 
            after finishing the job.
    """

    PROJECT_VERSION: ClassVar[str] = "20260728170528"
    jobname: str
    path: str
    results_npz_path: str
    pymoo_algorithm_bin_path: str
    initial_rng_state_path: str
    latest_rng_state: str
    prediction_pdbs_dir: str
    settings_json_path: str
    program_version_path: str
    predictor_input_dir: str
    predictor_output_dir: str
    profile_path: str
    essential_files: FrozenSet[str]



    @classmethod
    def create(cls, output_dir: str, jobname: Optional[str] = None) -> "PathsContainer":
        """
        Factory method to initialize the paths container for a given output directory.
        If no jobname is provided, it generates one based on the current timestamp.
        """
        abs_output_dir = os.path.abspath(output_dir)
        if jobname is None:
            today = datetime.today().strftime("%Y%m%d_%H%M%S")
            jobname = f"evodesign_{today}"
        base_path = os.path.join(abs_output_dir, jobname)
        results = os.path.join(base_path, "results.npz")
        pymoo_bin = os.path.join(base_path, "pymoo_algorithm.bin")
        init_rng = os.path.join(base_path, "initial_rng_state.txt")
        last_rng = os.path.join(base_path, "latest_rng_state.txt")
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
            latest_rng_state=last_rng,
            prediction_pdbs_dir=os.path.join(base_path, "prediction_pdbs"),
            settings_json_path=settings,
            program_version_path=version,
            predictor_input_dir=os.path.join(base_path, "predictor_input"),
            predictor_output_dir=os.path.join(base_path, "predictor_output"),
            profile_path=profile,
            essential_files=essential_files,
        )
