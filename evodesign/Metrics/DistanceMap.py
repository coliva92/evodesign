from .Metric import Metric
from .ContextInterface import ContextInterface
from Bio.PDB.Atom import Atom
import numpy as np
import numpy.typing as npt
from typing import List, Dict, Optional, Tuple
from .Normalization.Normalization import Normalization
from .Normalization.Reciprocal import Reciprocal
import itertools


class DistanceMap(Metric):

    def __init__(
        self,
        normalization: Optional[Normalization] = Reciprocal(),
    ):
        super().__init__()
        self.normalization = normalization

    @classmethod
    def compute_map(
        cls,
        ca_atoms: List[Atom],
    ) -> npt.NDArray[np.float64]:
        return np.array([a - b for (a, b) in itertools.combinations(ca_atoms, 2)])

    def do(
        self,
        model_map: npt.NDArray[np.float64],
        ref_map: npt.NDArray[np.float64],
        **kwargs,
    ) -> Tuple[float, Optional[float]]:
        # computing the mean over all individual values is the equivalent as computing
        # the weighted mean of the means of each value group
        rmse = np.sqrt(np.mean((ref_map - model_map) ** 2))
        norm = None
        if self.normalization is not None:
            norm = self.normalization.do(rmse)
        return rmse, norm

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> Dict[str, float]:
        model_ca_atoms = context.get_model_chain().ca_atoms
        ref_ca_atoms = context.get_reference_chain().ca_atoms
        model_map = self.compute_map(model_ca_atoms)
        ref_map = self.compute_map(ref_ca_atoms)
        rmse, norm = self.do(model_map, ref_map)
        return {
            "rmse": rmse,
            "norm_rmse": norm,
        }
