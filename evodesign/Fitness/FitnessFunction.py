from abc import ABC, abstractmethod
from typing import List
from ..RetrievableSettings import RetrievableSettings
import numpy.typing as npt
import numpy as np
from Bio.PDB.Structure import Structure
from Bio.PDB.Atom import Atom


class FitnessFunction(RetrievableSettings, ABC):

    def __init__(self, upper_bound: float) -> None:
        super().__init__()
        self.upper_bound = upper_bound

    def do(
        self,
        model_sequence: npt.NDArray[np.int64],
        model_structure: Structure,
        model_backbone: List[Atom],
        model_ca_atoms: List[Atom],
        ref_sequence: npt.NDArray[np.int64],
        ref_structure: Structure,
        ref_backbone: List[Atom],
        ref_ca_atoms: List[Atom],
    ) -> float:
        term_values = self.compute_term_values(
            model_sequence,
            model_structure,
            model_backbone,
            model_ca_atoms,
            ref_sequence,
            ref_structure,
            ref_backbone,
            ref_ca_atoms,
        )
        return self.compute_fitness(term_values)

    @abstractmethod
    def compute_term_values(
        self,
        model_sequence: npt.NDArray[np.int64],
        model_structure: Structure,
        model_backbone: List[Atom],
        model_ca_atoms: List[Atom],
        ref_sequence: npt.NDArray[np.int64],
        ref_structure: Structure,
        ref_backbone: List[Atom],
        ref_ca_atoms: List[Atom],
    ) -> npt.NDArray[np.float64]:
        raise NotImplementedError

    @abstractmethod
    def compute_fitness(self, term_values: npt.NDArray[np.float64]) -> float:
        raise NotImplementedError
