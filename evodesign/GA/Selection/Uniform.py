from .Selection import Selection
import numpy as np
import numpy.typing as npt


class Uniform(Selection):

    def do(
        self,
        rng: np.random.Generator,
        population: npt.NDArray[np.int64],
        fitness_values: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.int64]:
        parent_indices = np.concatenate(
            [
                rng.permutation(population.shape[0])
                for _ in range(self._num_parents_per_child)
            ]
        )
        return parent_indices
