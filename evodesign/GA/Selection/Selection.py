from abc import ABC, abstractmethod
from ...RetrievableSettings import RetrievableSettings
import numpy as np
import numpy.typing as npt


class Selection(RetrievableSettings, ABC):

    def __init__(self, num_parents_per_child: int = 2) -> None:
        super().__init__()
        self._num_parents_per_child = num_parents_per_child

    @abstractmethod
    def do(
        self,
        rng: np.random.Generator,
        population: npt.NDArray[np.int64],
        fitness_values: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.int64]:
        raise NotImplementedError
