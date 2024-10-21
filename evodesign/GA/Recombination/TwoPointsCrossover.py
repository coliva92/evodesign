from .Recombination import Recombination
import numpy as np
import numpy.typing as npt


class TwoPointsCrossover(Recombination):

    def _do(
        self,
        rng: np.random.Generator,
        mothers: npt.NDArray[np.int64],
        fathers: npt.NDArray[np.int64],
        offspring_mask: npt.NDArray[np._bool],
    ) -> npt.NDArray[np.int64]:
        crossover_points = rng.integers(0, mothers.shape[1], (mothers.shape[0], 2))
        sorting_mask = crossover_points[:, 0] > crossover_points[:, 1]
        crossover_points[sorting_mask] = crossover_points[sorting_mask][:, ::-1]
        offspring = mothers.copy()
        offspring[offspring_mask, crossover_points[0] : crossover_points[1] + 1] = (
            fathers[offspring_mask, crossover_points[0] : crossover_points[1] + 1]
        )
        return offspring
