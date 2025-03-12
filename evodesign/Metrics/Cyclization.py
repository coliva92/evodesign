from .Metric import Metric
from typing import List
from Bio.PDB.Atom import Atom
import evodesign.Utils.Normalization as Norm


class Cyclization(Metric):

    _mean = 1.3248119
    _deviation = 0.10498072

    def do(self, model_backbone: List[Atom], **kwargs) -> float:
        # assuming the backbone consists of atoms N-CA-C-O
        bond = model_backbone[-2] - model_backbone[0]
        return Norm.z_score(bond, self._mean, self._deviation)
