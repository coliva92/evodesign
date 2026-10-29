import os
from typing import List, Optional

from .AF3Interface import AF3Interface


class AlphaFold3Docker(AF3Interface):

    def __init__(
        self,
        docker_gpus: str = "all",  # can also admit MIG IDs; e.g.: 'device=MIG-blabla'
        using_rootless_docker: bool = True,
        jax_preallocate_memory: Optional[float] = 0.95,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.docker_gpus = docker_gpus
        self.using_rootless_docker = using_rootless_docker
        self.jax_preallocate_memory = jax_preallocate_memory
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

        cmd = ["docker", "run"]

        if not self.using_rootless_docker:
            uid = os.getuid()
            gid = os.getgid()
            cmd += [ "-u", f"{uid}:{gid}" ]

        if self.jax_preallocate_memory is None or self.jax_preallocate_memory <= 0:
            cmd += ["-e", "XLA_PYTHON_CLIENT_PREALLOCATE=false"]
        else:
            cmd += ["-e", f"XLA_CLIENT_MEM_FRACTION={self.jax_preallocate_memory}"]

        cmd += [
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
