from pymoo.core.survival import Survival
from pymoo.core.problem import Problem
from pymoo.core.population import Population
import numpy as np


class GenerationalElitism(Survival):

    def __init__(self, num_elites = 1) -> None:
        super().__init__(filter_infeasible=True)
        return

    def _do(
        self,
        problem: Problem,
        pop,
        n_survive=None,
        **kwargs,
    ):
        parents = pop[:n_survive]
        offspring = pop[n_survive:]
        parent_fitness = parents.get("F").flatten()
        best_idx = np.argmin(parent_fitness)
        # we pass [best_idx] so 'elite' remains a PyMOO Population instance
        elite = parents[[best_idx]] 
        next_generation = Population.merge(elite, offspring)
        return next_generation
