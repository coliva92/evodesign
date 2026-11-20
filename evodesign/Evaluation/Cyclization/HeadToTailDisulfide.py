from __future__ import annotations  # noqa: N999

from typing import ClassVar

import numpy as np
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue

from ...Protein.Geometry import find_atoms_in_residue
from ..ContextInterface import ContextInterface
from ..Normalization.Normalization import Normalization
from ..Normalization.Reciprocal import Reciprocal
from ..Normalization.ZScore import ZScore
from ..StructureBasedEvaluator import StructureBasedEvaluator


class HeadToTailDisulfide(StructureBasedEvaluator):
    # statistical values obtained from measuring relevant bond lengths
    # in a collection of cyclic peptide PDBs
    _Z_CALCULATORS: ClassVar[dict[str, Normalization]] = {
        "SG-SG": ZScore(
            mean_value=2.03,
            stdev=0.008,
        ),
        "SG-CB": ZScore(
            mean_value=3.023361386716186,
            stdev=0.04516770398282618,
        ),
    }
    _reciprocal: ClassVar[Reciprocal] = Reciprocal(scaling_factor=0.11111111111111111)

    def evaluate_bond_length(
        self,
        a: Atom,
        b: Atom,
        bond_name: str,
    ) -> tuple[float, float, float]:
        distance = a - b
        z_score = self._Z_CALCULATORS[bond_name].do(distance)
        norm_z_score = self._reciprocal(abs(z_score))
        return distance, z_score, norm_z_score

    def do(
        self,
        model_residues: list[Residue],
        **kwargs,
    ) -> tuple[float, ...]:
        amino_terminal = model_residues[0]
        tmp = find_atoms_in_residue(amino_terminal, ["CB", "SG"])
        cb1, s1 = tmp[0], tmp[1]
        carboxyl_terminal = model_residues[-1]
        tmp = find_atoms_in_residue(carboxyl_terminal, ["CB", "SG"])
        cb2, s2 = tmp[0], tmp[1]
        ss_distance, ss_z_score, ss_norm_z_score = self.evaluate_bond_length(
            s1,
            s2,
            "SG-SG",
        )
        h1_distance, h1_z_score, h1_norm_z_score = self.evaluate_bond_length(
            s1,
            cb2,
            "SG-CB",
        )
        h2_distance, h2_z_score, h2_norm_z_score = self.evaluate_bond_length(
            s2,
            cb1,
            "SG-CB",
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
    ) -> dict[str, float]:
        model_residues = context.get_model_chain().residues
        results = self.do(model_residues)
        return {
            "ss_cyclization": results[0],
            "ss_z_score": results[1],
            "ss_norm_z_score": results[2],
            "s1cb2_distance": results[3],
            "s1cb2_z_score": results[4],
            "s1cb2_norm_z_score": results[5],
            "s2cb1_distance": results[6],
            "s2cb1_z_score": results[7],
            "s2cb1_norm_z_score": results[8],
            "ss_cyclization_total_score": results[9],
        }
