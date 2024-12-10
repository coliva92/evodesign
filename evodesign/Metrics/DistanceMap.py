from .Metric import Metric
from typing import List
from Bio.PDB.Atom import Atom
import itertools
import numpy as np
import numpy.typing as npt


class DistanceMap(Metric):

    @classmethod
    def compute_map(cls, ca_atoms: List[Atom]) -> npt.NDArray[np.float64]:
        return np.array([a - b for (a, b) in itertools.combinations(ca_atoms, 2)])

    def do(
        self,
        model_map: npt.NDArray[np.float64],
        ref_map: npt.NDArray[np.float64],
        **kwargs
    ) -> float:
        # computing the mean over all individual values is the equivalent as computing
        # the weighted mean of the means of each value group
        rms = np.sqrt(np.mean((ref_map - model_map) ** 2))
        return rms
