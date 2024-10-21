from .Recombination import Recombination
import numpy as np
import numpy.typing as npt


class SinglePointCrossover(Recombination):

    def _do(
        self,
        rng: np.random.Generator,
        mothers: npt.NDArray[np.int64],
        fathers: npt.NDArray[np.int64],
        offspring_mask: npt.NDArray[np._bool],
    ) -> npt.NDArray[np.int64]:
        crossover_point = rng.integers(0, mothers.shape[1] + 1, mothers.shape[0])
        offspring = mothers.copy()
        offspring[offspring_mask, crossover_point:] = fathers[
            offspring_mask, crossover_point:
        ]
        return offspring
