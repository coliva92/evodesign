from .StructuralMetric import StructuralMetric
from .ContextInterface import ContextInterface
from Bio.PDB.Atom import Atom
from typing import List, Dict, Tuple
from scipy.spatial import cKDTree
import numpy as np
import numpy.typing as npt


class lDDT(StructuralMetric):

    def __init__(
        self,
        ca_atoms_only: bool = False,
        radius: float = 15,
        cutoffs: List[float] = [0.5, 1.0, 2.0, 4.0],
    ) -> None:
        super().__init__()
        self.ca_atoms_only = ca_atoms_only
        self.radius = radius
        self.cutoffs = cutoffs
        self._num_atoms_per_residue = 4 if not self.ca_atoms_only else 1
        return

    def do(
        self,
        model_distances: Dict[Tuple[int, int], float],
        ref_distances: Dict[Tuple[int, int], float],
        **kwargs,
    ) -> npt.NDArray[np.float64]:
        per_residue_lddt = []
        for i in range(len(ref_distances)):
            lddt = self.compute_lddt(model_distances[i], ref_distances[i])
            per_residue_lddt.append(lddt)
        return np.array(per_residue_lddt)

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> Dict[str, float]:
        ref_distances_dict = context.get_extra_param_value("ref_distances_dict")
        if ref_distances_dict is None:
            ref_atoms = (
                context.get_reference_chain().ca_atoms
                if self.ca_atoms_only
                else context.get_reference_chain().backbone_atoms
            )
            ref_distances_dict = self.get_distances_dict(ref_atoms)
            context.set_extra_param_value("ref_distances_dict", ref_distances_dict)
        model_atoms = (
            context.get_model_chain().ca_atoms
            if self.ca_atoms_only
            else context.get_model_chain().backbone_atoms
        )
        model_distances_dict = self.get_distances_dict(model_atoms)
        per_residue_lddt = self.do(model_distances_dict, ref_distances_dict)
        return {"lddt": per_residue_lddt.mean()}

    def get_distances_dict(
        self,
        atoms: List[Atom],
    ) -> Dict[Tuple[int, int], float]:
        coordinates = self._get_coordinates_list(atoms)
        tree = cKDTree(coordinates)
        neighborhoods = tree.query_ball_point(coordinates, self.radius)
        distances = (len(neighborhoods) // self._num_atoms_per_residue) * [None]
        for center_idx, neighbors in enumerate(neighborhoods):
            residue_idx = center_idx // self._num_atoms_per_residue
            distances[residue_idx] = {}
            center_atom = atoms[center_idx]
            # remove neighboring atoms belonging to the same residue as center
            for i in neighbors:
                neighbor_atom = atoms[i]
                if self._belong_to_same_residue(neighbor_atom, center_atom):
                    continue
                distances[residue_idx][(center_idx, i)] = center_atom - neighbor_atom
        return distances

    def compute_lddt(
        self,
        model_distances: Dict[Tuple[int, int], float],
        ref_distances: Dict[Tuple[int, int], float],
    ) -> float:
        # remove from the model distances those that do not appear in the reference
        working_distances = {}
        for atoms_pair in model_distances:
            if atoms_pair in ref_distances:
                working_distances[atoms_pair] = model_distances[atoms_pair]
        L = len(ref_distances)
        ratios_sum = 0
        for c in self.cutoffs:
            M = 0
            for atoms_pair, model_dist in working_distances.items():
                ref_dist = ref_distances[atoms_pair]
                if abs(ref_dist - model_dist) <= c:
                    M += 1
            ratios_sum += M / L
        return ratios_sum / len(self.cutoffs)

    def _get_coordinates_list(
        self,
        atoms: List[Atom],
    ) -> npt.NDArray[np.float32]:
        return np.array([atom.get_coord() for atom in atoms])

    def _belong_to_same_residue(
        self,
        a: Atom,
        b: Atom,
    ) -> bool:
        return a.get_full_id()[3][1] == b.get_full_id()[3][1]
