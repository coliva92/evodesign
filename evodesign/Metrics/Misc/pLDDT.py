import numpy as np  # noqa: N999

from ...Chemistry.Chain import Chain
from ..ContextInterface import ContextInterface
from ..StructureBasedEvaluator import StructureBasedEvaluator


class pLDDT(StructureBasedEvaluator):
    def do(
        self,
        model_chain: Chain,
        **kwargs,
    ) -> tuple[float, ...]:
        bfactors = tuple(
            atom.get_bfactor()
            if atom.get_bfactor() <= 1.0
            else atom.get_bfactor() / 100.0
            for atom in model_chain.structure.get_atoms()
        )
        return bfactors

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> dict[str, float]:
        model_chain = context.get_model_chain()
        per_residue_plddt = self.do(model_chain)
        return {"plddt": np.mean(per_residue_plddt)}
