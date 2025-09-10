from abc import ABC, abstractmethod
from .Predictor import Predictor
from .DirectoryManager import DirectoryManager


class ESMFoldInterface(Predictor, ABC):

    def predict_single_pdb_file(
        self,
        sequence: str,
        directory: DirectoryManager,
    ) -> None:
        prediction = self.predict_single_pdb_str(sequence)
        with open(directory.prediction_pdbs_dir, "wt", encoding="utf-8") as pdb_file:
            pdb_file.write(prediction)

    @abstractmethod
    def predict_single_pdb_str(
        self,
        sequence: str,
    ) -> str:
        raise NotImplementedError
