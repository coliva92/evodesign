from pymoo.core.problem import Problem
from pymoo.core.survival import Survival


class Generational(Survival):

    def __init__(self):
        super().__init__(filter_infeasible=False)

    def _do(
        self,
        problem: Problem,
        pop,
        n_survive=None,
        **kwargs,
    ):
        return pop[n_survive:]
