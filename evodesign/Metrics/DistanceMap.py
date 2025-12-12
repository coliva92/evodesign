from .StructuralMetric import StructuralMetric
from .ContextInterface import ContextInterface
from Bio.PDB.Atom import Atom
import numpy as np
import numpy.typing as npt
from typing import List, Dict, Optional, Tuple
from .Normalization.Normalization import Normalization
from .Normalization.Reciprocal import Reciprocal
import itertools


class DistanceMap(StructuralMetric):

    def __init__(
        self,
        ca_atoms_only: bool = True,
        normalization: Optional[Normalization] = Reciprocal(),
    ):
        super().__init__()
        self.ca_atoms_only = ca_atoms_only
        self.normalization = normalization
        return

    @classmethod
    def compute_map(
        cls,
        atoms: List[Atom],
    ) -> npt.NDArray[np.float64]:
        return np.array([a - b for (a, b) in itertools.combinations(atoms, 2)])

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
        ref_map = context.get_extra_param_value("ref_distance_map")
        if ref_map is None:
            ref_atoms = (
                context.get_reference_chain().ca_atoms
                if self.ca_atoms_only
                else context.get_reference_chain().backbone_atoms
            )
            ref_map = self.compute_map(ref_atoms)
            context.set_extra_param_value("ref_distance_map", ref_map)
        model_atoms = (
            context.get_model_chain().ca_atoms
            if self.ca_atoms_only
            else context.get_model_chain().backbone_atoms
        )
        model_map = self.compute_map(model_atoms)
        rmse, norm = self.do(model_map, ref_map)
        return {
            "rmse": rmse,
            "norm_rmse": norm,
        }
