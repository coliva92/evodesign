from pymoo.core.survival import Survival
from pymoo.core.problem import Problem


class Generational(Survival):

    def __init__(self) -> None:
        super().__init__(filter_infeasible=False)

    def _do(self, problem: Problem, pop, n_survive=None, **kwargs):
        return pop[n_survive:]
