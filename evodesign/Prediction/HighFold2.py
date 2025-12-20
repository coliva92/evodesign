from .AlphaFoldInterface import AlphaFoldInterface
from typing import Optional, List
import os


class HighFold2(AlphaFoldInterface):

    def __init__(
        self,
        path_to_prediction_py: str,
        stop_at_score: Optional[float] = None,
        num_recycle: int = 0,
        num_models: int = 1,
        host_url: Optional[str] = None,
        data_dir: Optional[str] = None,
        msa_mode: str = "single_sequence",  # {mmseqs2_uniref_env,mmseqs2_uniref,single_sequence}
        use_dropout: bool = True,
        use_gpu_relax: bool = False,
        # disulfide_bond_pairs: Optional[list] = None,
        # index_ss: Optional[list] = None,
        flag_cyclic_peptide: bool = True,
        flag_nc: bool = False,
        # custom_dist_cst: Optional[list] = None,
    ) -> None:
        super().__init__()
        self.path_to_prediction_py = os.path.abspath(path_to_prediction_py)
        self.stop_at_score = stop_at_score
        self.num_recycle = num_recycle
        self.num_models = num_models
        self.host_url = host_url
        self.data_dir = os.path.abspath(data_dir)
        self.msa_mode = msa_mode
        self.use_dropout = use_dropout
        self.use_gpu_relax = use_gpu_relax
        # self.disulfide_bond_pairs = disulfide_bond_pairs
        # self.index_ss = index_ss
        self.flag_cyclic_peptide = flag_cyclic_peptide
        self.flag_nc = flag_nc
        # self.custom_dist_cst = custom_dist_cst

    def _create_model_input(
        self,
        sequence: str,
        protein_full_name: str,
        input_dir: str,
        output_dir: str,
    ) -> str:
        fasta_path = os.path.join(input_dir, f"{protein_full_name}.fasta")
        with open(fasta_path, "wt", encoding="utf-8") as fasta_file:
            fasta_file.write(f">{protein_full_name}\n{sequence}\n")
        return fasta_path

    def _prediction_pdb_path(
        self,
        protein_full_name: str,
        output_dir: str,
    ) -> str:
        return os.path.join(output_dir, protein_full_name, "ranked_0.pdb")

    def _create_cmd_array(
        self,
        input_path: str,
        output_dir: str,
        do_batch_inference: bool,
    ) -> List[str]:
        cmd = [
            "python3",
            self.path_to_prediction_py,
            input_path,
            output_dir,
            f"--num-recycle {self.num_recycle}",
            f"--num-models {self.num_models}",
            f"--msa-mode {self.msa_mode}",
        ]
        if self.stop_at_score is not None:
            cmd.append(f"--stop-at-score {self.stop_at_score}")
        if self.host_url is not None:
            cmd.append(f"--host-url {self.stop_at_score}")
        if self.data_dir is not None:
            cmd.append(f"--data {self.data_dir}")
        if self.use_dropout:
            cmd.append(f"--use-dropout")
        if self.use_gpu_relax:
            cmd.append(f"--use_gpu_relax")
        if self.flag_cyclic_peptide:
            cmd.append(f"--flag-cyclic-peptide")
        if self.flag_nc:
            cmd.append(f"--flag-nc")
        return cmd
