import os
from abc import ABC, abstractmethod

from .DirectoryManager import DirectoryManager
from .Predictor import Predictor




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
        return



    @abstractmethod
    def predict_single_pdb_str(
        self,
        sequence: str,
    ) -> str:
        raise NotImplementedError
