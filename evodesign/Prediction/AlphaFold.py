from .Predictor import Predictor
from ..Utils.Subprocess import run_subprocess
import os
import shutil


class AlphaFold(Predictor):

    def __init__(
        self,
        path_to_create_fakemsa_py: str,
        path_to_run_alphafold_py: str,
        output_dir: str,
        mgnify_database_path: str,
        data_dir: str,
        max_template_date: str = "2020-05-14",
        model_preset: str = "monomer",
        db_preset: str = "reduced_dbs",
    ) -> None:
        super().__init__()
        self.path_to_create_fakemsa_py = os.path.abspath(path_to_create_fakemsa_py)  # https://github.com/Zuricho/ParaFold_dev/blob/main/parafold/create_fakemsa.py
        self.path_to_run_alphafold_py = os.path.abspath(path_to_run_alphafold_py)
        self.mgnify_database_path = os.path.abspath(mgnify_database_path)  # /media/biocomp/My\ Passport/mgnify/mgy_clusters_2018_12.fa
        self.data_dir = os.path.abspath(data_dir)  # /media/biocomp/My\ Passport/reduced_dbs
        self.max_template_date = max_template_date
        self.model_preset = (
            model_preset  # { 'monomer', 'monomer_casp14', 'monomer_ptm', 'multimer' }
        )
        self.db_preset = db_preset  # { 'reduced_dbs', 'full_dbs' }
        self.output_dir = os.path.dirname(os.path.abspath(output_dir))

    def predict_pdb_str(self, sequence: str) -> str:
        pdb_path = os.path.join(self.output_dir, "prediction.pdb.tmp")
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
        # run the script for creating an empty MSA
        run_subprocess(
            [
                "python3",
                self.path_to_create_fakemsa_py,
                f"--fasta_paths={fasta_path}",
                f"--output_dir={self.output_dir}",
            ]
        )
        self.run_alphafold_docker(fasta_path)
        prediction_pdb = os.path.join(self.output_dir, protein_name, "ranked_0.pdb")
        shutil.copyfile(prediction_pdb, pdb_path)
        os.remove(fasta_path)
        self.delete_output_dir()

    def run_alphafold_docker(self, fasta_path: str):
        run_subprocess(
            [
                "python3",
                self.path_to_run_alphafold_py,
                "--use_precomputed_msas=True",
                f"--fasta_paths={fasta_path}",
                f"--max_template_date={self.max_template_date}",
                f"--model_preset={self.model_preset}",
                f"--db_preset={self.db_preset}",
                f"--output_dir={self.output_dir}",
                f"--mgnify_database_path={self.mgnify_database_path}",
                f"--data_dir={self.data_dir}",
            ]
        )

    def delete_output_dir(self) -> None:
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)
