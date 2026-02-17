from .StructuralMetric import StructuralMetric
from .ContextInterface import ContextInterface
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from typing import List, Tuple, Dict
import evodesign.Utils.Normalization as Norm
from ..Utils.Geometry import find_atoms_in_residue, compute_dihedral_angle
import numpy as np


class LengthDisulfideCyclization(StructuralMetric):

    _bond_lengths = {
        "SG-SG": 2.03,
        "SG-CB": 3.018093966326573,       
    }
    _bond_stdevs = {
        "SG-SG": 0.008,
        "SG-CB": 0.07019316291606081,
    }
    _scaling_factors = {
        "SG-SG": 0.11111111111111111,
        "SG-CB": 0.11111111111111111,
    }

    def evaluate_bond_length(
        self, a: Atom, b: Atom, bond_name: str
    ) -> Tuple[float, float, float]:
        distance = a - b
        avg = self._bond_lengths[bond_name]
        stdev = self._bond_stdevs[bond_name]
        z_score = Norm.z_score(distance, avg, stdev)
        norm_z_score = Norm.reciprocal(abs(z_score), self._scaling_factors[bond_name])
        return (distance, z_score, norm_z_score)

    def do(
        self,
        model_residues: List[Residue],
        **kwargs,
    ) -> Tuple[float]:
        amino_terminal = model_residues[0]
        tmp = find_atoms_in_residue(amino_terminal, ["CB", "SG"])
        cb1, s1 = tmp[0], tmp[1]
        carboxyl_terminal = model_residues[-1]
        tmp = find_atoms_in_residue(carboxyl_terminal, ["CB", "SG"])
        cb2, s2 = tmp[0], tmp[1]
        ss_distance, ss_z_score, ss_norm_z_score = self.evaluate_bond_length(
            s1, s2, "SG-SG"
        )
        h1_distance, h1_z_score, h1_norm_z_score = self.evaluate_bond_length(
            s1, cb2, "SG-CB"
        )
        h2_distance, h2_z_score, h2_norm_z_score = self.evaluate_bond_length(
            s2, cb1, "SG-CB"
        )
        ss_total_score = np.mean(
            [
                ss_norm_z_score,
                h1_norm_z_score,
                h2_norm_z_score,
            ]
        )
        return (
            ss_distance,
            ss_z_score,
            ss_norm_z_score,
            h1_distance,
            h1_z_score,
            h1_norm_z_score,
            h2_distance,
            h2_z_score,
            h2_norm_z_score,
            ss_total_score,
        )

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> Dict[str, float]:
        model_residues = context.get_model_chain().residues
        results = self.do(model_residues)
        return {
            "ss_cyclization": results[0],
            "ss_z_score": results[1],
            "ss_norm_z_score": results[2],
            "h1_distance": results[3],
            "h1_z_score": results[4],
            "h1_norm_z_score": results[5],
            "h2_distance": results[6],
            "h2_z_score": results[7],
            "h2_norm_z_score": results[8],
            "ss_cyclization_total_score": results[9],
        }
