import os
from typing import List

from .AF3Interface import AF3Interface


class AlphaFold3Docker(AF3Interface):

    def __init__(
        self,
        docker_gpus: str = "all",  # can also admit MIG IDs; e.g.: 'device=MIG-blabla'
        **kwargs,
    ):
        self.docker_gpus = docker_gpus
        return

    def _create_cmd_array(
        self,
        input_path: str,
        output_dir: str,
        do_batch_inference: bool,
    ) -> List[str]:
        input_path = os.path.realpath(input_path)
        output_dir = os.path.realpath(output_dir)
        model_dir = os.path.realpath(self.model_dir)

        uid = os.getuid()
        gid = os.getgid()

        cmd = [
            "docker",
            "run",
            "-u",
            f"{uid}:{gid}",
            "--volume",
            (
                f"{input_path}:/data/af_input"
                if do_batch_inference
                else f"{input_path}:/data/af_input"
            ),
            "--volume",
            f"{output_dir}:/data/af_output",
            "--volume",
            f"{model_dir}:/data/models",
            "--gpus",
            self.docker_gpus,
            "alphafold3",
            "python",
            "run_alphafold.py",
            (
                f"--input_dir=/data/af_input"
                if do_batch_inference
                else f"--json_path=/data/af_input/input.json"
            ),
            "--model_dir=/data/models",
            "--output_dir=/data/af_output",
        ]
        cmd.extend(self._get_af3_flags())
        return cmd
