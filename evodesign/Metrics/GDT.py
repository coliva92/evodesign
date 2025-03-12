from .Metric import Metric
from typing import List, Optional
import numpy as np
from Bio.PDB.Atom import Atom
from Bio.PDB.Superimposer import Superimposer


class GDT(Metric):

    def __init__(self, cutoffs: List[float] = [1.0, 2.0, 4.0, 8.0]) -> None:
        self.cutoffs = cutoffs

    def do(
        self,
        model_backbone: List[Atom],
        ref_backbone: List[Atom],
        superimposer: Optional[Superimposer] = None,
        **kwargs
    ) -> float:
        # assume the backbones are superimposed if no superimposer is provided
        if superimposer is not None:
            superimposer.set_atoms(ref_backbone, model_backbone)
            superimposer.apply(model_backbone)
        distances = np.array([a - b for a, b in zip(model_backbone, ref_backbone)])
        gdt = np.mean([np.mean([d <= c for d in distances]) for c in self.cutoffs])
        return gdt
