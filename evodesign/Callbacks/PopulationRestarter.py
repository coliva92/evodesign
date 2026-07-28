import numpy as np
from pymoo.algorithms.base.genetic import GeneticAlgorithm
from pymoo.core.algorithm import Algorithm as PyMOOAlgorithm
from pymoo.core.callback import Callback
from pymoo.core.population import Population

from ..Statistics import (ALPHABET_SIZE, get_population_amino_acid_loss,
                          get_population_identity)


class PopulationRestarter(Callback):

    def __init__(
        self,
        diversity_loss_tol: float = 0.95,
    ) -> None:
        super().__init__()
        self.diversity_loss_tol = diversity_loss_tol
        return

    def notify(self, algorithm: PyMOOAlgorithm) -> None:
        assert isinstance(algorithm, GeneticAlgorithm)
        population = np.array([solution for solution in algorithm.pop.get("X")])
        d1 = get_population_amino_acid_loss(population) / ALPHABET_SIZE
        d2 = get_population_identity(population) / population.shape[1]
        similarity = (d1 + d2) / 2
        if similarity >= self.diversity_loss_tol:
            # get the elite individual
            fitness = algorithm.pop.get("F").flatten()
            best_idx = np.argmin(fitness)
            elite = algorithm.pop[[best_idx]]
            # generate new random sequences
            n_pop = len(algorithm.pop)
            n_new = n_pop - 1
            new_pop = algorithm.initialization.sampling(
                algorithm.problem, n_new, algorithm=algorithm
            )
            # we must evaluate the individuals before returning to the algorithm
            algorithm.evaluator.eval(algorithm.problem, new_pop, algorithm=algorithm)
            algorithm.pop = Population.merge(elite, new_pop)
        return
