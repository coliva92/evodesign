import itertools  # noqa: N999

import numpy as np
import numpy.typing as npt
from Bio.PDB.Atom import Atom

from ..ContextInterface import ContextInterface
from ..Normalization.Normalization import Normalization
from ..StructuralMetric import StructuralMetric


class DistanceMap(StructuralMetric):
    def __init__(
        self,
        ca_atoms_only: bool = True,
        normalization: Normalization | None = None,
    ):
        super().__init__()
        self.ca_atoms_only = ca_atoms_only
        self.normalization = normalization

    @classmethod
    def compute_pairwise_distances(
        cls,
        atoms: list[Atom],
    ) -> npt.NDArray[np.float64]:
        return np.array([a - b for (a, b) in itertools.combinations(atoms, 2)])

    def do(
        self,
        model_distances: npt.NDArray[np.float64],
        ref_distances: npt.NDArray[np.float64],
        **kwargs,
    ) -> tuple[float, ...]:
        # computing the mean over all individual values is the equivalent as computing
        # the weighted mean of the means of each value group
        rmse = np.sqrt(np.mean((ref_distances - model_distances) ** 2))
        norm = None
        if self.normalization is not None:
            norm = self.normalization.do(rmse)
        return rmse, norm

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> dict[str, float]:
        ref_distances = context.get_extra_param_value("ref_distance_map")
        if ref_distances is None:
            ref_atoms = (
                context.get_reference_chain().ca_atoms
                if self.ca_atoms_only
                else context.get_reference_chain().backbone_atoms
            )
            ref_distances = self.compute_map(ref_atoms)
            context.set_extra_param_value("ref_distance_map", ref_distances)
        model_atoms = (
            context.get_model_chain().ca_atoms
            if self.ca_atoms_only
            else context.get_model_chain().backbone_atoms
        )
        model_distances = self.compute_pairwise_distances(model_atoms)
        rmse, norm = self.do(model_distances, ref_distances)
        return {
            "rmse": rmse,
            "norm_rmse": norm,
        }
