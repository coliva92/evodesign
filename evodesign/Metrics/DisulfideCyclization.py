from .StructuralMetric import StructuralMetric
from .ContextInterface import ContextInterface
from Bio.PDB.Residue import Residue
from typing import List, Tuple, Dict
import evodesign.Utils.Normalization as Norm


class DisulfideCyclization(StructuralMetric):

    _mean = 2.04 # from Chaney MO & Steinrauf LK (1974) Acta Cryst B30:711-716
    _stdev = 0.006 # from natural disulfides
    _scaling_factor = 0.05

    def do(
        self,
        model_residues: List[Residue],
        **kwargs,
    ) -> Tuple[float]:
        amino_terminal = model_residues[0]
        carboxyl_terminal = model_residues[-1]
        a, b = None, None
        for atom in amino_terminal.get_atoms():
            name = atom.get_name()
            if name == "SG" or name == "S":
                a = atom
                break
        for atom in carboxyl_terminal.get_atoms():
            name = atom.get_name()
            if name == "SG" or name == "S":
                b = atom
                break
        distance = a - b
        z_score = Norm.z_score(distance, self._mean, self._stdev)
        norm_z_score = Norm.reciprocal(z_score, self._scaling_factor)
        return (distance, z_score, norm_z_score)

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> Dict[str, float]:
        model_residues = context.get_model_chain().residues
        distance, z_score, norm_z_score = self.do(model_residues)
        return {
            "ss_cyclization": distance,
            "z_score": z_score,
            "norm_z_score": norm_z_score,
        }
