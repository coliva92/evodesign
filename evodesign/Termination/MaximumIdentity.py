from pymoo.termination.max_gen import MaximumGenerationTermination
from ..Utils.Statistics import get_population_amino_acid_loss, get_population_identity
import numpy as np


class MaximumIdentity(MaximumGenerationTermination):

    def __init__(
        self,
        n_max_gen: float = float("inf"),
        max_similarity: float = 1.0,
    ) -> None:
        super().__init__(n_max_gen)
        self.n_max_gen = n_max_gen
        self.max_similarity = max_similarity
        return

    def _update(self, algorithm):
        progress = super()._update(algorithm)
        if progress >= 1.0:
            return progress
        population = np.array([
            solution
            for solution in algorithm.pop.get("X")
        ])
        d1 = get_population_amino_acid_loss(population)
        d2 = get_population_identity(population, sample_size=30)
        similarity = (d1 + d2) / 2
        return 1.0 if similarity >= self.max_similarity else progress
