from .Metric import Metric
import numpy as np
from typing import Optional, List
from Bio.PDB.Atom import Atom
from Bio.PDB.Superimposer import Superimposer
from .ContextInterface import ContextInterface


class TMScore(Metric):

    def normalizing_constant(self, n: int) -> float:
        return 1.24 * np.cbrt(n - 15) - 1.8

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
        d0 = self.normalizing_constant(len(ref_backbone))
        distances = np.array([a - b for a, b in zip(model_backbone, ref_backbone)])
        tm_score = np.mean(1 / (1 + (distances / d0) ** 2))
        return tm_score

    def do_for_fitness_fn(self, context: ContextInterface):
        model_backbone = context.get_model_chain().backbone_atoms
        ref_backbone = context.get_model_chain().backbone_atoms
        superimposer = context.get_extra_param_value("superimposer")
        if superimposer is None:
            superimposer = Superimposer()
            context.set_extra_param_value("superimposer", superimposer)
        tm_score = self.do(model_backbone, ref_backbone, superimposer)
        return {"tm_score": tm_score}
