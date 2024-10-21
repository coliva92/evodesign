from .Selection import Selection
import numpy as np
import numpy.typing as npt


class Tournament(Selection):

    def __init__(self, tournament_size: int) -> None:
        super().__init__()
        self.tournament_size = tournament_size

    def do(
        self,
        rng: np.random.Generator,
        population: npt.NDArray[np.int64],
        fitness_values: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.int64]:
        num_permutations = self._num_parents_per_child * self.tournament_size
        permutations = np.concatenate(
            [rng.permutation(population.shape[0]) for _ in range(num_permutations)]
        )
        tournament_groups = permutations.reshape(
            self._num_parents_per_child * population.shape[0], self.tournament_size
        )
        parent_indices = np.apply_along_axis(
            lambda group_indices: group_indices[fitness_values[group_indices].argmax()],
            1,
            tournament_groups,
        )
        return parent_indices
