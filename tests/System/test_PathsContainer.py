import os
from evodesign.System.PathsContainer import PathsContainer




def test_paths_container_creates_expected_paths(tmp_path):
    """
    Test that PathsContainer correctly resolves all internal paths 
    when given an output directory and a specific jobname.
    """
    jobname = "test_job"
    container = PathsContainer.create(output_dir=str(tmp_path), jobname=jobname)
    
    base_path = os.path.join(str(tmp_path), jobname)
    
    assert container.jobname == jobname
    assert container.path == base_path
    assert container.results_npz_path == os.path.join(base_path, "results.npz")
    assert container.pymoo_algorithm_bin_path == os.path.join(base_path, "pymoo_algorithm.bin")
    assert container.initial_rng_state_path == os.path.join(base_path, "initial_rng_state.txt")
    assert container.latest_rng_state == os.path.join(base_path, "latest_rng_state.txt")
    assert container.settings_json_path == os.path.join(base_path, "settings.json")
    assert container.program_version_path == os.path.join(base_path, "version.txt")
    assert container.profile_path == os.path.join(base_path, "aa_profile.txt")
    assert container.prediction_pdbs_dir == os.path.join(base_path, "prediction_pdbs")
    assert container.predictor_input_dir == os.path.join(base_path, "predictor_input")
    assert container.predictor_output_dir == os.path.join(base_path, "predictor_output")
    
    # Check essential files
    expected_essential_files = {
        container.results_npz_path,
        container.pymoo_algorithm_bin_path,
        container.initial_rng_state_path,
        container.latest_rng_state,
        container.settings_json_path,
        container.program_version_path,
        container.profile_path
    }
    assert container.essential_files == expected_essential_files


def test_paths_container_generates_timestamp_jobname_if_none(tmp_path):
    """
    Test that PathsContainer generates a jobname with the current timestamp 
    when no jobname is explicitly provided.
    """
    container = PathsContainer.create(output_dir=str(tmp_path))
    
    assert container.jobname is not None
    assert container.jobname.startswith("evodesign_")
    
    # Verify that the generated jobname is used consistently for path resolution
    base_path = os.path.join(str(tmp_path), container.jobname)
    assert container.path == base_path
    assert container.results_npz_path == os.path.join(base_path, "results.npz")
