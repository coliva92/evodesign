from .Metric import Metric
from ..Utils.Chain import Chain
from .ContextInterface import ContextInterface
from typing import Dict
import numpy as np


class pLDDT(Metric):

    def uses_predictor(self) -> bool:
        return True

    def do(
        self,
        model_chain: Chain,
        **kwargs,
    ) -> float:
        bfactors = np.array(
            [
                (
                    atom.get_bfactor()
                    if atom.get_bfactor() <= 1.0
                    else atom.get_bfactor() / 100.0
                )
                for atom in model_chain.structure.get_atoms()
            ]
        )
        return bfactors.mean()

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> Dict[str, float]:
        model_chain = context.get_model_chain()
        plddt = self.do(model_chain)
        return {"plddt": plddt}
