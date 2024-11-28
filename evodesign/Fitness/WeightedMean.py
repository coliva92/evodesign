from .FitnessFunction import FitnessFunction
import numpy as np
from typing import List
import numpy.typing as npt
from Bio.PDB.Structure import Structure
from Bio.PDB.Atom import Atom


class WeightedMean(FitnessFunction):

    _terms = {
        "rmsd_weight": 0,
        "gdt_weight": 1,
        "tm_score_weight": 2,
        "distance_map_weight": 3,
        "cyclization_weight": 4,
        "rosetta_energy_weight": 5,
        "plddt_weight": 6,
    }

    def __init__(
        self,
        upper_bound: float = 1.0,
        rmsd_weight: float = 1.0,
        gdt_weight: float = 1.0,
        tm_score_weight: float = 1.0,
        distance_map_weight: float = 1.0,
        cyclization_weight: float = 1.0,
        rosetta_energy_weight: float = 1.0,
        plddt_weight: float = 1.0,
    ) -> None:
        super().__init__(upper_bound)
        self.rmsd_weight = rmsd_weight
        self.gdt_weight = gdt_weight
        self.tm_score_weight = tm_score_weight
        self.distance_map_weight = distance_map_weight
        self.cyclization_weight = cyclization_weight
        self.rosetta_energy_weight = rosetta_energy_weight
        self.plddt_weight = plddt_weight
        self._weights = np.array(
            [
                self.rmsd_weight,
                self.gdt_weight,
                self.tm_score_weight,
                self.distance_map_weight,
                self.cyclization_weight,
                self.rosetta_energy_weight,
                self.plddt_weight,
            ]
        )

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
        # TODO después de refactorizar los metrics, incluirlos aquí
        return np.array(
            [rmsd, gdt, tm_score, distance_map_rms, cyc_z_score, energy_norm, plddt]
        )

    def compute_fitness(self, term_values: npt.NDArray[np.float64]) -> float:
        return np.average(term_values, weights=self._weights)

    def set_weight(self, term_name: str, weight: float) -> None:
        self._weights[self._terms[term_name]] = self.__dict__[term_name] = weight
