from .Metric import Metric
from ..Utils.Chain import Chain
from .ContextInterface import ContextInterface
from typing import Dict


class pLDDT(Metric):

    def do(self, model_chain: Chain, **kwargs) -> float:
        plddt = 0  # calculate from model_chain
        return plddt

    def do_for_fitness_fn(self, context: ContextInterface) -> Dict[str, float]:
        model_chain = context.get_model_chain()
        plddt = self.do(model_chain)
        return {"plddt": plddt}
