from .Predictor import Predictor
from typing import List, Tuple
from Bio.PDB.Atom import Atom


class Null(Predictor):

    def do(self, sequence: str, pdb_path: str) -> Tuple[List[Atom], float]:
        return [], 0.0

    def predict_pdb_str(self, sequence: str) -> str:
        pass
