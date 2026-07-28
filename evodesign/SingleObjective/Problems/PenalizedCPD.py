import numpy as np

from ...Chemistry.Sequences import compute_identity
from .SingleObjectiveCPD import SingleObjectiveCPD




class PenalizedCPD(SingleObjectiveCPD):

    def __init__(self, identity_tol: float = 0.3, **kwargs):
        super().__init__(n_ieq_constr=1, **kwargs)
        self.identity_tol = identity_tol
        return



    def _evaluate_constraints(
        self,
        x,
        out,
        *args,
        **kwargs,
    ) -> None:
        population_size = x.shape[0]
        out["G"] = np.array(
            [
                compute_identity(x[i, :], self.ref_chain.sequence_numpy)
                - self.identity_tol
                for i in range(population_size)
            ]
        )
        return



    def _evaluate(
        self,
        x,
        out,
        *args,
        **kwargs,
    ) -> None:
        super()._evaluate(x, out, *args, **kwargs)
        self._evaluate_constraints(x, out, *args, **kwargs)
        return
