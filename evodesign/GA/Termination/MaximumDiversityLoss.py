from pymoo.termination.max_gen import MaximumGenerationTermination
from pymoo.core.algorithm import Algorithm as PyMOOAlgorithm
from ...Statistics import get_population_diversity_loss
from typing import Optional
import numpy as np


class MaximumDiversityLoss(MaximumGenerationTermination):

    def __init__(
        self,
        n_max_gen: float = float("inf"),
        diversity_loss_tol: float = 1.0,
        sample_size: Optional[int] = None,
    ) -> None:
        super().__init__(n_max_gen)
        self.n_max_gen = n_max_gen
        self.diversity_loss_tol = diversity_loss_tol
        self.sample_size = sample_size
        return

    def _update(self, algorithm: PyMOOAlgorithm):
        progress = super()._update(algorithm)
        if progress >= 1.0 or algorithm.n_iter % 10 != 0:
            return progress
        population = np.array([solution for solution in algorithm.pop.get("X")])
        diversity_loss = get_population_diversity_loss(population, self.sample_size)
        return 1.0 if diversity_loss >= self.diversity_loss_tol else progress
