from pymoo.termination.max_gen import MaximumGenerationTermination
from ...Statistics import (
    get_population_amino_acid_loss,
    get_population_identity,
    ALPHABET_SIZE,
)
from typing import Optional
import numpy as np


class MaximumDiversityLoss(MaximumGenerationTermination):

    def __init__(
        self,
        n_max_gen: float = float("inf"),
        max_similarity: float = 1.0,
        sample_size: Optional[int] = None,
    ) -> None:
        super().__init__(n_max_gen)
        self.n_max_gen = n_max_gen
        self.max_similarity = max_similarity
        self.sample_size = sample_size
        return

    def _update(self, algorithm):
        progress = super()._update(algorithm)
        if progress >= 1.0 or algorithm.n_iter % 10 != 0:
            return progress
        population = np.array([solution for solution in algorithm.pop.get("X")])
        d1 = get_population_amino_acid_loss(population) / ALPHABET_SIZE
        d2 = get_population_identity(population, self.sample_size) / population.shape[1]
        similarity = (d1 + d2) / 2
        return 1.0 if similarity >= self.max_similarity else progress
