from typing import List

from ...System.Subprocess import run_subprocess
from .AlphaFold3Base import AlphaFold3Base




class AlphaFold3Docker(AlphaFold3Base):

    def _create_cmd_array(
        self,
        input_path: str,
        output_dir: str,
        do_batch_inference: bool,
    ) -> List[str]:
        cmd = [
            "docker",
            "run",
            "-it",
            "--volume",
            (
                f"{input_path}:/root/af_input"
                if do_batch_inference
                else f"{input_path}:/root/af_input/input.json"
            ),
            "--volume",
            f"{output_dir}:/root/af_output",
            "--volume",
            f"{self.model_dir}:/root/models",
            "--gpus",
            "all",
            "alphafold3",
            "python",
            "run_alphafold.py",
            (
                f"--input_dir=/root/af_input"
                if do_batch_inference
                else f"--json_path=/root/af_input/input.json"
            ),
            "--model_dir=/root/models",
            "--output_dir=/root/af_output",
        ]
        cmd.extend(self._get_af3_flags())
        return cmd



    def run_inference(
        self,
        input_path: str,
        output_dir: str,
        do_batch_inference: bool,
    ) -> None:
        run_subprocess(
            self._create_cmd_array(input_path, output_dir, do_batch_inference)
        )
        # AF3's docker image creates the output files as root;
        # we must change the owner
        cmd = [
            "docker",
            "run",
            "-it",
            "--volume",
            f"{output_dir}:/root/af_output",
            "alphafold3",
            "chown",
            "-R",
            "1000:1000",
            "/root/af_output",
        ]
        run_subprocess(cmd)
        return
