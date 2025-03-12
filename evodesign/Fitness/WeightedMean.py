from .FitnessFunction import FitnessFunction
import numpy as np
from typing import List
import numpy.typing as npt
from Bio.PDB.Structure import Structure
from Bio.PDB.Atom import Atom
from ..Metrics.RMSD import RMSD
from ..Metrics.GDT import GDT
from ..Metrics.TMScore import TMScore
from ..Metrics.DistanceMap import DistanceMap
from ..Metrics.Cyclization import Cyclization
from ..Metrics.RosettaEnergyFunction import RosettaEnergyFunction


class WeightedMean(FitnessFunction):

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
        rmsd_calc: RMSD = RMSD(),
        gdt_calc: GDT = GDT(),
        tm_score_calc: TMScore = TMScore(),
        dm_calc: DistanceMap = DistanceMap(),
        cyc_calc: Cyclization = Cyclization(),
        eng_calc: RosettaEnergyFunction = RosettaEnergyFunction(),
    ) -> None:
        super().__init__(upper_bound)
        self.rmsd_weight = rmsd_weight
        self.gdt_weight = gdt_weight
        self.tm_score_weight = tm_score_weight
        self.distance_map_weight = distance_map_weight
        self.cyclization_weight = cyclization_weight
        self.rosetta_energy_weight = rosetta_energy_weight
        self.plddt_weight = plddt_weight
        self.rmsd_calc = rmsd_calc
        self.gdt_calc = gdt_calc
        self.tm_score_calc = tm_score_calc
        self.dm_calc = dm_calc
        self.cyc_calc = cyc_calc
        self.eng_calc = eng_calc
        self._weights = np.array(
            [weight for name, weight in self.__dict__.items() if "_weight" in name]
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
        kwargs = locals()
        return np.array(
            [
                calc.do(**kwargs)
                for name, calc in self.__dict__.items()
                if "_calc" in name
            ]
        )

    def compute_fitness(self, term_values: npt.NDArray[np.float64]) -> float:
        return np.average(term_values, weights=self._weights)
