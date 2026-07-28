import pytest

from evodesign.Prediction.AlphaFold3 import AlphaFold3
from evodesign.Prediction.AlphaFold3Docker import AlphaFold3Docker
from evodesign.Prediction.HighFold3 import HighFold3
from evodesign.Prediction.HighFold3Docker import HighFold3Docker


def test_alphafold3_docker_cmd_array():
    predictor = AlphaFold3Docker(
        model_dir="/path/to/models",
        num_recycles=3,
        num_diffusion_samples=5,
        version=2,
        run_data_pipeline=True,
    )
    cmd = predictor._create_cmd_array("/input", "/output", do_batch_inference=True)
    
    assert cmd[0] == "docker"
    assert cmd[1] == "run"
    assert "--volume" in cmd
    assert "alphafold3" in cmd
    assert "python" in cmd
    assert "run_alphafold.py" in cmd
    assert "--run_data_pipeline=True" in cmd
    assert "--num_recycles=3" in cmd
    assert "--num_diffusion_samples=5" in cmd
    assert "--force_output_dir=True" in cmd


def test_alphafold3_local_cmd_array():
    predictor = AlphaFold3(
        path_to_run_alphafold_py="/path/to/run_alphafold.py",
        model_dir="/path/to/models",
        num_recycles=3,
        num_diffusion_samples=5,
        version=2,
        run_data_pipeline=True,
    )
    cmd = predictor._create_cmd_array("/input", "/output", do_batch_inference=True)
    
    assert cmd[0] == "python3"
    assert cmd[1] == "/path/to/run_alphafold.py"
    assert "--input_dir=/input" in cmd
    assert "--model_dir=/path/to/models" in cmd
    assert "--output_dir=/output" in cmd
    assert "--run_data_pipeline=True" in cmd
    assert "--num_recycles=3" in cmd
    assert "--num_diffusion_samples=5" in cmd
    assert "--force_output_dir=True" in cmd
    assert "docker" not in cmd


def test_highfold3_docker_cmd_array():
    predictor = HighFold3Docker(
        model_dir="/path/to/models",
        head_to_tail=False,
        num_recycles=2,
        run_data_pipeline=False,
    )
    cmd = predictor._create_cmd_array("/input", "/output", do_batch_inference=True)
    
    assert cmd[0] == "docker"
    assert cmd[1] == "run"
    assert "highfold3" in cmd
    assert "python" in cmd
    assert "run_alphafold.py" in cmd
    assert "--run_data_pipeline=False" in cmd
    assert "--num_recycles=2" in cmd
    assert "--head_to_tail=False" in cmd
    # HighFold3 should remove --force_output_dir
    assert "--force_output_dir=True" not in cmd


def test_highfold3_local_cmd_array():
    predictor = HighFold3(
        path_to_run_alphafold_py="/path/to/run_alphafold.py",
        model_dir="/path/to/models",
        head_to_tail=True,
        disulfide_chain_res=[[1, 2], [3, 4]],
        num_recycles=2,
        run_data_pipeline=False,
    )
    cmd = predictor._create_cmd_array("/input", "/output", do_batch_inference=True)
    
    assert cmd[0] == "python3"
    assert cmd[1] == "/path/to/run_alphafold.py"
    assert "--run_data_pipeline=False" in cmd
    assert "--num_recycles=2" in cmd
    assert "--head_to_tail=True" in cmd
    assert "--disulfide_chain_res [[1, 2], [3, 4]]" in cmd
    assert "--force_output_dir=True" not in cmd
    assert "docker" not in cmd
