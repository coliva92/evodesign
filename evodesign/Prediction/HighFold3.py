from .AlphaFold3 import AlphaFold3
from typing import List, Optional


class HighFold3(AlphaFold3):

    def __init__(
        self,
        path_to_run_alphafold_py: str,
        model_dir: str,
        output_dir: str,
        num_recycles: int = 1,
        num_diffusion_samples: int = 1,
        run_data_pipeline: bool = False,
        model_seeds: Optional[List[int]] = None,
        version: int = 2,
        head_to_tail: bool = True,
        disulfide_chain_res: Optional[List[List[int]]] = None,
    ):
        super().__init__(
            path_to_run_alphafold_py,
            model_dir,
            output_dir,
            num_recycles,
            num_diffusion_samples,
            run_data_pipeline,
            model_seeds,
            version,
        )
        self.head_to_tail = head_to_tail
        self.disulfide_chain_res = disulfide_chain_res

    def _get_cmd_array(self, json_path: str) -> List[str]:
        cmd = super()._get_cmd_array(json_path)
        cmd.append(f"--head_to_tail={self.head_to_tail}")
        if self.disulfide_chain_res is not None:
            cmd.append(f"--disulfide_chain_res {self.disulfide_chain_res}")
        return cmd
