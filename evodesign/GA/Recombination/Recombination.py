from abc import ABC, abstractmethod
from ...RetrievableSettings import RetrievableSettings
import numpy as np
import numpy.typing as npt


class Recombination(RetrievableSettings, ABC):

    def __init__(
        self,
        recombination_probability: float = 1.0,
        num_parents_per_offspring: int = 2,
        num_offspring_per_parents: int = 1,
    ) -> None:
        super().__init__()
        self.recombination_probability = recombination_probability
        self._num_parents_per_offspring = num_parents_per_offspring
        self._num_offspring_per_parents = num_offspring_per_parents

    def do(
        self,
        rng: np.random.Generator,
        population: npt.NDArray[np.int64],
        parent_indices: npt.NDArray[np.int64],
    ) -> npt.NDArray[np.int64]:
        if parent_indices.shape[0] % self._num_parents_per_offspring != 0:
            raise RuntimeError
        num_offspring = parent_indices.shape[0] // self._num_parents_per_offspring
        recombine_indices = rng.random(num_offspring) < self.recombination_probability
        mothers = population[parent_indices[::2]]
        fathers = population[parent_indices[1::2]]
        return self._do(rng, mothers, fathers, recombine_indices)

    @abstractmethod
    def _do(
        self,
        rng: np.random.Generator,
        mothers: npt.NDArray[np.int64],
        fathers: npt.NDArray[np.int64],
        recombine_indices: npt.NDArray[np.bool_]
    ) -> npt.NDArray[np.int64]:
        raise NotImplementedError
