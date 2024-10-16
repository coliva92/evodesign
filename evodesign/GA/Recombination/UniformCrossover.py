from .Recombination import Recombination
from typing import List
import numpy as np
import numpy.typing as npt


class UniformCrossover(Recombination):

    def _do(
        self,
        rng: np.random.Generator,
        mothers: npt.NDArray[np.int64],
        fathers: npt.NDArray[np.int64],
        recombine_indices: npt.NDArray[np.bool_],
    ) -> npt.NDArray[np.int64]:
        offspring = mothers.copy()
        mask = rng.integers(0, 2, mothers.shape)
        complement = 1 - mask
        offspring[recombine_indices] = (
            mask[recombine_indices] * mothers[recombine_indices]
            + complement[recombine_indices] * fathers[recombine_indices]
        )
        return offspring
