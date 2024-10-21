from .Recombination import Recombination
import numpy as np
import numpy.typing as npt


class MiddlePointCrossover(Recombination):

    def _do(
        self,
        rng: np.random.Generator,
        mothers: npt.NDArray[np.int64],
        fathers: npt.NDArray[np.int64],
        offspring_mask: npt.NDArray[np._bool],
    ) -> npt.NDArray[np.int64]:
        crossover_point = mothers.shape[1] // 2
        offspring = mothers.copy()
        offspring[offspring_mask, crossover_point:] = fathers[
            offspring_mask, crossover_point:
        ]
        return offspring
