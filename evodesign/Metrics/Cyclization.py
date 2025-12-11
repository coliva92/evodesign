from .Metric import Metric
from .ContextInterface import ContextInterface
from typing import List, Tuple, Dict
from Bio.PDB.Atom import Atom
import evodesign.Utils.Normalization as Norm


class Cyclization(Metric):

    _mean = 1.3248119
    _stdev = 0.10498072
    _scaling_factor = 0.05

    def uses_predictor(self) -> bool:
        return True

    def do(
        self,
        model_backbone: List[Atom],
        **kwargs,
    ) -> Tuple[float]:
        # assuming the backbone consists of atoms N-CA-C-O
        distance = model_backbone[-2] - model_backbone[0]
        z_score = Norm.z_score(distance, self._mean, self._stdev)
        norm_z_score = Norm.reciprocal(z_score, self._scaling_factor)
        return (distance, z_score, norm_z_score)

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> Dict[str, float]:
        model_backbone = context.get_model_chain().backbone_atoms
        distance, z_score, norm_z_score = self.do(model_backbone)
        return {
            "cyclization": distance,
            "z_score": z_score,
            "norm_z_score": norm_z_score,
        }
