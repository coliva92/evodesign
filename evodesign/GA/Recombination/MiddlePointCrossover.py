from .Recombination import Recombination
import numpy as np
import numpy.typing as npt


class MiddlePointCrossover(Recombination):

    def _do(
        self,
        rng: np.random.Generator,
        mothers: npt.NDArray[np.int64],
        fathers: npt.NDArray[np.int64],
        recombine_indices: npt.NDArray[np._bool],
    ) -> npt.NDArray[np.int64]:
        offspring = mothers.copy()
        crossover_point = mothers.shape[1] // 2
        temp_chromosomes = np.hstack(
            (mothers[:crossover_point], fathers[crossover_point:])
        )
        offspring[recombine_indices] = temp_chromosomes[recombine_indices]
        return offspring
