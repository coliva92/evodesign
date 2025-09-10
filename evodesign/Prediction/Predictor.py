from abc import ABC, abstractmethod
from ..RetrievableSettings import RetrievableSettings
from .DirectoryManager import DirectoryManager
from typing import List
import os
from ..Utils.Exceptions import HttpGatewayTimeout
from requests.exceptions import ConnectTimeout


class Predictor(RetrievableSettings, ABC):

    @abstractmethod
    def predict_single_pdb_file(
        self,
        sequence: str,
        directory: DirectoryManager,
    ) -> None:
        raise NotImplementedError

    @abstractmethod
    def do(
        self,
        sequences: List[str],
        directory: DirectoryManager,
    ) -> None:
        # The default behavior is to predict the PDB files of all the sequences
        # in the population.
        # Also, PDBParser cannot produce an instance of `Structure` directly from
        # a raw PDB string, it can only do it by reading from a PDB file.
        # Thus, the predicted structure must be stored first, before
        # loading the `Structure` instance.
        os.makedirs(directory.prediction_pdbs_dir, exist_ok=True)
        directory.empty_folders_content()
        for i, sequence in enumerate(sequences):
            protein_name = f"{directory.prefix}_{i}"
            pdb_path = os.path.join(
                directory.prediction_pdbs_dir, f"{protein_name}.pdb"
            )
            while True:
                try:
                    self.predict_single_pdb_file(sequence, pdb_path, directory)
                    break
                except (HttpGatewayTimeout, ConnectTimeout):
                    continue
