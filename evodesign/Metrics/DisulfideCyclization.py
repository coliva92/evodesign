from .StructuralMetric import StructuralMetric
from .ContextInterface import ContextInterface
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from typing import List, Tuple, Dict
import evodesign.Utils.Normalization as Norm
from ..Utils.Geometry import find_atoms_in_residue, compute_dihedral_angle
import numpy as np


class DisulfideCyclization(StructuralMetric):

    _bond_lengths = {
        "SG-SG": 2.0,  # from Oliviero C (2022) Life 12(7):986
        "CB-CB": 3.8,  # from Srinivasan N et al, (1990) Int J Peptide Protein Res 36:147-155
    }
    _bond_stdevs = {
        "SG-SG": 0.1,
        "CB-CB": 0.2,
    }
    _torsion_angles = {
        "X1": np.deg2rad(60),  # from Deplazes E et al. (2019) Proteins 88(3):485-502
        "X3": np.deg2rad(90),  # from Deplazes E et al. (2019) Proteins 88(3):485-502
    }
    _torsion_stdevs = {
        "X1": np.deg2rad(30),
        "X3": np.deg2rad(30),
    }
    _scaling_factors = {"SG-SG": 0.05, "CB-CB": 0.05, "X1": 0.05, "X3": 0.05}

    def evaluate_bond_length(self, a: Atom, b: Atom, bond_name: str) -> Tuple[float]:
        distance = a - b
        avg = self._bond_lengths[bond_name]
        stdev = self._bond_stdevs[bond_name]
        z_score = Norm.z_score(distance, avg, stdev)
        norm_z_score = Norm.reciprocal(abs(z_score), self._scaling_factors[bond_name])
        return (distance, z_score, norm_z_score)

    def evaluate_dihedral_angle(
        self,
        a: Atom,
        b: Atom,
        c: Atom,
        d: Atom,
        torsion_name: str,
    ) -> Tuple[float]:
        dihedral = compute_dihedral_angle(a, b, c, d)
        avg = self._torsion_angles[torsion_name]
        stdev = self._torsion_stdevs[torsion_name]
        z_score = Norm.z_score(abs(dihedral), avg, stdev)
        norm_z_score = Norm.reciprocal(
            abs(z_score), self._scaling_factors[torsion_name]
        )
        return (dihedral, z_score, norm_z_score)

    def do(
        self,
        model_residues: List[Residue],
        **kwargs,
    ) -> Tuple[float]:
        amino_terminal = model_residues[0]
        tmp = find_atoms_in_residue(amino_terminal, ["C", "CA", "CB", "SG"])
        c1, ca1, cb1, s1 = tmp[0], tmp[1], tmp[2], tmp[3]
        carboxyl_terminal = model_residues[-1]
        tmp = find_atoms_in_residue(carboxyl_terminal, ["C", "CA", "CB", "SG"])
        c2, ca2, cb2, s2 = tmp[0], tmp[1], tmp[2], tmp[3]
        ss_distance, ss_z_score, ss_norm_z_score = self.evaluate_bond_length(
            s1, s2, "SG-SG"
        )
        cb_distance, cb_z_score, cb_norm_z_score = self.evaluate_bond_length(
            cb1, cb2, "CB-CB"
        )
        x1_a, x1_a_z_score, x1_a_norm_z_score = self.evaluate_dihedral_angle(
            c1, ca1, cb1, s1, "X1"
        )
        x1_b, x1_b_z_score, x1_b_norm_z_score = self.evaluate_dihedral_angle(
            c2, ca2, cb2, s2, "X1"
        )
        x3, x3_z_score, x3_norm_z_score = self.evaluate_dihedral_angle(
            cb1, s1, s2, cb2, "X3"
        )
        total_score = np.mean(
            [
                ss_norm_z_score,
                cb_norm_z_score,
                x1_a_norm_z_score,
                x1_b_norm_z_score,
                x3_norm_z_score,
            ]
        )
        return (
            ss_distance,
            ss_z_score,
            ss_norm_z_score,
            cb_distance,
            cb_z_score,
            cb_norm_z_score,
            np.rad2deg(x1_a),
            x1_a_z_score,
            x1_a_norm_z_score,
            np.rad2deg(x1_b),
            x1_b_z_score,
            x1_b_norm_z_score,
            np.rad2deg(x3),
            x3_z_score,
            x3_norm_z_score,
            total_score,
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
            "cb_bond": results[3],
            "cb_z_score": results[4],
            "cb_norm_z_score": results[5],
            "torsion_a": results[6],
            "torsion_a_z_score": results[7],
            "torsion_a_norm_z_score": results[8],
            "torsion_b": results[9],
            "torsion_b_z_score": results[10],
            "torsion_b_norm_z_score": results[11],
            "torsion_ss": results[12],
            "torsion_ss_z_score": results[13],
            "torsion_ss_norm_z_score": results[14],
            "ss_cyclization_total_score": results[15]
        }
