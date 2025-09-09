from abc import ABC, abstractmethod
from typing import List, Tuple
from ..Utils.Chain import ChainFactory, Chain
from ..RetrievableSettings import RetrievableSettings
import numpy as np
import os
import shutil
from ..Utils.Exceptions import HttpGatewayTimeout
from requests.exceptions import ConnectTimeout


class Predictor(RetrievableSettings, ABC):

    @abstractmethod
    def predict_single_pdb_str(self, sequence: str) -> str:
        raise NotImplementedError

    def predict_single_pdb_file(self, sequence: str, pdb_path: str) -> None:
        if not os.path.isfile(pdb_path):
            os.makedirs(os.path.dirname(os.path.abspath(pdb_path)), exist_ok=True)
        prediction = self.predict_single_pdb_str(sequence)
        with open(pdb_path, "wt", encoding="utf-8") as pdb_file:
            pdb_file.write(prediction)

    def do(
        self, sequences: List[str], pdbs_dir: str, pdb_prefix: str = "tmp_prediction"
    ) -> None:
        # The default behavior is to predict the PDB files of all the sequences
        # in the population.
        # Also, PDBParser cannot produce an instance of `Structure` directly from
        # a raw PDB string, it can only do it by reading from a PDB file.
        # Thus, the predicted structure must be stored first, before
        # loading the `Structure` instance.
        os.makedirs(pdbs_dir, exist_ok=True)
        for i, sequence in enumerate(sequences):
            pdb_path = os.path.join(pdbs_dir, f"{pdb_prefix}_{i}.pdb")
            while True:
                try:
                    self.predict_single_pdb_file(sequence, pdb_path)
                    break
                except (HttpGatewayTimeout, ConnectTimeout):
                    continue

    def compute_plddt(self, chain: Chain) -> float:
        bfactors = np.array(
            [
                (
                    atom.get_bfactor()
                    if atom.get_bfactor() <= 1.0
                    else atom.get_bfactor() / 100.0
                )
                for atom in chain.structure.get_atoms()
            ]
        )
        return bfactors.mean()

    def delete_folder(self, folder_dir: str) -> None:
        if os.path.exists(folder_dir):
            shutil.rmtree(folder_dir)
