from abc import ABC, abstractmethod
from ...RetrievableSettings import RetrievableSettings
import numpy as np
import numpy.typing as npt
from typing import Optional


class Mutation(RetrievableSettings, ABC):

    def __init__(self, mutation_probability: float = 1.0) -> None:
        super().__init__()
        self.mutation_probability = mutation_probability

    @abstractmethod
    def do(
        self,
        rng: np.random.Generator,
        children: npt.NDArray[np.int64],
        constraints: Optional[npt.NDArray[np.float64]] = None,
    ) -> npt.NDArray[np.int64]:
        raise NotImplementedError
