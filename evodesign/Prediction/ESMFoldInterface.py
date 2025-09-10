from abc import ABC, abstractmethod
from .Predictor import Predictor
from .DirectoryManager import DirectoryManager
import os


class ESMFoldInterface(Predictor, ABC):

    def predict_single_pdb_file(
        self,
        sequence: str,
        protein_name: str,
        directory: DirectoryManager,
    ) -> None:
        prediction = self.predict_single_pdb_str(sequence)
        pdb_path = os.path.join(
            directory.prediction_pdbs_dir, f"{directory.prefix}_{protein_name}.pdb"
        )
        with open(pdb_path, "wt", encoding="utf-8") as pdb_file:
            pdb_file.write(prediction)

    @abstractmethod
    def predict_single_pdb_str(
        self,
        sequence: str,
    ) -> str:
        raise NotImplementedError
