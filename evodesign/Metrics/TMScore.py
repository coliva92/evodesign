from .StructuralMetric import StructuralMetric
import numpy as np
from typing import Optional, List, Dict, Tuple
from Bio.PDB.Atom import Atom
from Bio.PDB.Superimposer import Superimposer
from .ContextInterface import ContextInterface
from .Normalization.Normalization import Normalization


class TMScore(StructuralMetric):

    def __init__(
        self,
        residue_weights: Optional[List[float]] = None,
        normalization: Optional[Normalization] = None,
    ) -> None:
        super().__init__()
        self.residue_weights = residue_weights
        self.normalization = normalization
        self._backbone_weights = None
        return

    def normalizing_constant(
        self,
        n: int,
    ) -> float:
        return 1.24 * np.cbrt(n - 15) - 1.8

    def do(
        self,
        model_backbone: List[Atom],
        ref_backbone: List[Atom],
        superimposer: Optional[Superimposer] = None,
        **kwargs,
    ) -> Tuple[float, float]:
        if self.residue_weights is not None and self._backbone_weights is None:
            self._backbone_weights = []
            for i in range(len(self.residue_weights)):
                w = self.residue_weights[i] / 4.0
                self._backbone_weights.extend(4 * [w])
        # assume the backbones are superimposed if no superimposer is provided
        if superimposer is not None:
            superimposer.set_atoms(ref_backbone, model_backbone)
            superimposer.apply(model_backbone)
        d0 = self.normalizing_constant(len(ref_backbone))
        distances = np.array([a - b for a, b in zip(model_backbone, ref_backbone)])
        tmp = 1 / (1 + (distances / d0) ** 2)
        tm_score = np.average(tmp, weights=self._backbone_weights)
        norm_tm_score = tm_score
        if self.normalization is not None:
            norm_tm_score = self.normalization.do(tm_score)
        return tm_score, norm_tm_score

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> Dict[str, float]:
        model_backbone = context.get_model_chain().backbone_atoms
        ref_backbone = context.get_reference_chain().backbone_atoms
        superimposer = context.get_extra_param_value("superimposer")
        if superimposer is None:
            superimposer = Superimposer()
            context.set_extra_param_value("superimposer", superimposer)
        tm_score, norm_tm_score = self.do(model_backbone, ref_backbone, superimposer)
        return {"tm_score": tm_score, "norm_tm_score": norm_tm_score}
