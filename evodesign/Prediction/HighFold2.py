from .Predictor import Predictor
from ..Utils.Subprocess import run_subprocess
from typing import Optional
import os
import shutil


class HighFold2(Predictor):

    def __init__(
        self,
        path_to_prediction_py: str,
        output_dir: str,
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
        self.output_dir = os.path.dirname(os.path.abspath(output_dir))
        self.stop_at_score = stop_at_score
        self.num_recycle = num_recycle
        self.num_models = num_models
        self.host_url = host_url
        self.data_dir = data_dir
        self.msa_mode = msa_mode
        self.use_dropout = use_dropout
        self.use_gpu_relax = use_gpu_relax
        # self.disulfide_bond_pairs = disulfide_bond_pairs
        # self.index_ss = index_ss
        self.flag_cyclic_peptide = flag_cyclic_peptide
        self.flag_nc = flag_nc
        # self.custom_dist_cst = custom_dist_cst

    def predict_pdb_str(self, sequence: str) -> str:
        pdb_path = "tmp_prediction.pdb"
        self.predict_pdb_file(sequence, pdb_path)
        with open(pdb_path, "rt", encoding="utf-8") as pdb_file:
            prediction = pdb_file.read()
        os.remove(pdb_path)
        return prediction

    def predict_pdb_file(self, sequence: str, pdb_path: str) -> None:
        os.makedirs(self.output_dir, exist_ok=True)
        protein_name = os.path.splitext(os.path.basename(pdb_path))[0]
        fasta_path = os.path.join(self.output_dir, f"{protein_name}.fasta")
        with open(fasta_path, "wt", encoding="utf-8") as fasta_file:
            fasta_file.write(f">{protein_name}\n{sequence}\n")
        self.run_highfold(fasta_path)
        prediction_pdb = os.path.join(self.output_dir, protein_name, "ranked_0.pdb")
        shutil.copyfile(prediction_pdb, pdb_path)
        os.remove(fasta_path)
        self.delete_output_dir()

    def run_highfold(self, fasta_path: str):
        cmd = [
            "python3",
            self.path_to_prediction_py,
            fasta_path,
            self.output_dir,
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
        run_subprocess(cmd)

    def delete_output_dir(self) -> None:
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)
