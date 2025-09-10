from .Predictor import Predictor
from .DirectoryManager import DirectoryManager
from typing import List


class Null(Predictor):

    def predict_single_pdb_file(
        self,
        sequence: str,
        directory: DirectoryManager,
    ) -> None:
        pass

    def do(
        self,
        sequences: List[str],
        directory: DirectoryManager,
    ) -> None:
        return
