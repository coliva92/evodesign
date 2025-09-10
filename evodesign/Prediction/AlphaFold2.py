from .AlphaFoldInterface import AlphaFoldInterface
from ..Utils.Subprocess import run_subprocess
import os


class AlphaFold2(AlphaFoldInterface):

    def __init__(
        self,
        path_to_create_fakemsa_py: str,
        path_to_run_alphafold_py: str,
        mgnify_database_path: str,
        data_dir: str,
        max_template_date: str = "2020-05-14",
        model_preset: str = "monomer",
        db_preset: str = "reduced_dbs",
    ):
        super().__init__()
        self.path_to_create_fakemsa_py = os.path.abspath(
            path_to_create_fakemsa_py
        )  # https://github.com/Zuricho/ParaFold_dev/blob/main/parafold/create_fakemsa.py
        self.path_to_run_alphafold_py = os.path.abspath(path_to_run_alphafold_py)
        self.mgnify_database_path = os.path.abspath(
            mgnify_database_path
        )  # /media/biocomp/My\ Passport/mgnify/mgy_clusters_2018_12.fa
        self.data_dir = os.path.abspath(
            data_dir
        )  # /media/biocomp/My\ Passport/reduced_dbs
        self.max_template_date = max_template_date
        self.model_preset = (
            model_preset  # { 'monomer', 'monomer_casp14', 'monomer_ptm', 'multimer' }
        )
        self.db_preset = db_preset  # { 'reduced_dbs', 'full_dbs' }

    def _create_model_input(
        self,
        sequence: str,
        protein_name: str,
        input_dir: str,
        output_dir: str,
    ) -> str:
        fasta_path = os.path.join(input_dir, f"{protein_name}.fasta")
        with open(fasta_path, "wt", encoding="utf-8") as fasta_file:
            fasta_file.write(f">{protein_name}\n{sequence}\n")
        self.create_empty_msa(fasta_path, output_dir)
        return fasta_path

    def _prediction_pdb_path(
        self,
        output_dir: str,
        protein_name: str,
    ) -> str:
        return os.path.join(output_dir, protein_name, "ranked_0.pdb")

    def _create_cmd_array(
        self,
        input_path: str,
        output_dir: str,
        do_batch_inference: bool,
    ):
        return [
            "python3",
            self.path_to_run_alphafold_py,
            "--use_precomputed_msas=True",
            f"--fasta_paths={input_path}",
            f"--max_template_date={self.max_template_date}",
            f"--model_preset={self.model_preset}",
            f"--db_preset={self.db_preset}",
            f"--output_dir={output_dir}",
            f"--mgnify_database_path={self.mgnify_database_path}",
            f"--data_dir={self.data_dir}",
        ]

    def create_empty_msa(
        self,
        fasta_path: str,
        output_dir: str,
    ) -> None:
        # run the script for creating an empty MSA
        run_subprocess(
            [
                "python3",
                self.path_to_create_fakemsa_py,
                f"--fasta_paths={fasta_path}",
                f"--output_dir={output_dir}",
            ]
        )
