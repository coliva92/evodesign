from pymoo.core.survival import Survival
from pymoo.core.problem import Problem
from pymoo.core.population import Population
import numpy as np


class GenerationalElitism(Survival):

    def __init__(self) -> None:
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
        parent_fitness = parents.get("F")
        best_idx = np.argmin(parent_fitness[:, 0])
        # we pass [best_idx] so 'elite' remains a PyMOO Population instance
        elite = parents[[best_idx]] 
        offspring_fitness = offspring.get("F")
        worst_idx = np.argmax(offspring_fitness[:, 0])
        # create a boolean mask to remove the worst individual
        keep_mask = np.ones(len(offspring), dtype=bool)
        keep_mask[worst_idx] = False
        surviving_offspring = offspring[keep_mask]
        next_generation = Population.merge(elite, surviving_offspring)
        return next_generation[:n_survive]
