from .Predictor import Predictor
from ..Utils.Subprocess import run_subprocess
from Bio.PDB import MMCIFParser, PDBIO
import os
import shutil
import json
from typing import List, Optional


class AlphaFold3(Predictor):

    _parser = MMCIFParser()
    _io = PDBIO()

    def __init__(
        self,
        path_to_run_alphafold_py: str,
        model_dir: str,
        output_dir: str,
        # db_dir: str
        run_data_pipeline: bool = False,
        model_seeds: Optional[List[int]] = None,
        version: int = 3,
    ):
        super().__init__()
        self.path_to_run_alphafold_py = os.path.abspath(path_to_run_alphafold_py)
        self.model_dir = os.path.abspath(model_dir)
        self.output_dir = os.path.abspath(output_dir)
        # self.db_dir = db_dir
        self.run_data_pipeline = run_data_pipeline
        if model_seeds is None:
            import time

            model_seeds = [int(time.time())]
        self.model_seeds = model_seeds
        self.version = version

    def _get_cmd_array(self, json_path: str) -> List[str]:
        return [
            "python3",
            self.path_to_run_alphafold_py,
            f"--json_path={json_path}",
            f"--model_dir={self.model_dir}",
            f"--output_dir={self.output_dir}",
            f"--run_data_pipeline={self.run_data_pipeline}",
            "--force_output_dir=True",
        ]

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
        input_json = {
            "name": protein_name,
            "modelSeeds": self.model_seeds,
            "sequences": [
                {
                    "protein": {
                        "id": "A",
                        "sequence": sequence,
                        "unpairedMsa": "",
                        "pairedMsa": "",
                        "templates": [],
                    }
                }
            ],
            "dialect": "alphafold3",
            "version": self.version,
        }
        json_path = os.path.join(self.output_dir, f"{protein_name}.json")
        with open(json_path, "wt", encoding="utf-8") as json_file:
            json.dump(input_json, json_file)
        self.run_alphafold(json_path)
        prediction_cif = os.path.join(
            self.output_dir, protein_name, f"{protein_name}_model.cif"
        )
        prediction_pdb = self.convert_cif_to_pdb(prediction_cif, "A")
        shutil.copyfile(prediction_pdb, pdb_path)
        os.remove(json_path)
        self.delete_output_dir()

    def run_alphafold(self, json_path: str):
        run_subprocess(self._get_cmd_array(json_path))

    def convert_cif_to_pdb(self, cif_path: str, structure_id: str) -> str:
        structure = self._parser.get_structure(structure_id, cif_path)
        pdb_path = cif_path.replace(".cif", ".pdb")
        self._io.set_structure(structure)
        self._io.save(pdb_path)
        return pdb_path

    def delete_output_dir(self) -> None:
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)
