from abc import ABC, abstractmethod
from ...RetrievableSettings import RetrievableSettings
import numpy as np
import numpy.typing as npt


class Recombination(RetrievableSettings, ABC):

    def __init__(
        self,
        recombination_probability: float = 1.0,
        num_parents_per_child: int = 2,
        num_offspring_per_parents_group: int = 1,
    ) -> None:
        super().__init__()
        self.recombination_probability = recombination_probability
        self._num_parents_per_child = num_parents_per_child
        self._num_offspring_per_parents_group = num_offspring_per_parents_group
    
    @abstractmethod
    def get_offspring(
        self,
        rng: np.random.Generator,
        mothers: npt.NDArray[np.int64],
        fathers: npt.NDArray[np.int64],
        offspring_mask: npt.NDArray[np.bool_],
    ) -> npt.NDArray[np.int64]:
        raise NotImplementedError

    def do(
        self,
        rng: np.random.Generator,
        population: npt.NDArray[np.int64],
        parent_indices: npt.NDArray[np.int64],
    ) -> npt.NDArray[np.int64]:
        if parent_indices.shape[0] % self._num_parents_per_child != 0:
            raise RuntimeError
        num_offspring = parent_indices.shape[0] // self._num_parents_per_child
        offspring_mask = rng.random(num_offspring) < self.recombination_probability
        mothers = population[parent_indices[::2]]
        fathers = population[parent_indices[1::2]]
        return self.get_offspring(rng, mothers, fathers, offspring_mask)
