from .Predictor import Predictor
from typing import List, Tuple
from ..Utils.Chain import Chain


class Null(Predictor):

    def do(
        self, sequences: List[str], pdbs_dir: str, pdb_prefix: str = "tmp_prediction"
    ) -> None:
        return

    def predict_single_pdb_str(self, sequence: str) -> str:
        pass

    def predict_single_pdb_file(self, sequence: str, pdb_path: str) -> None:
        pass
