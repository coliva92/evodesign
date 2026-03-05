from pymoo.core.callback import Callback
from pymoo.core.algorithm import Algorithm as PyMOOAlgorithm
from pymoo.core.population import Population
from ..Statistics import (
    get_population_amino_acid_loss,
    get_population_identity,
    ALPHABET_SIZE,
)
from typing import Optional
import numpy as np


class PopulationRestarter(Callback):

    def __init__(
        self,
        diversity_loss_tol: float = 0.95,
        next_callback: Optional[Callback] = None,
    ) -> None:
        super().__init__()
        self.next_callback = next_callback
        self.diversity_loss_tol = diversity_loss_tol
        return

    def notify(self, algorithm: PyMOOAlgorithm) -> None:
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
        if self.next_callback is not None:
            self.next_callback.notify(algorithm)
        return
