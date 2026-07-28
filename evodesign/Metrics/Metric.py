from abc import ABC, abstractmethod
from typing import Dict, Tuple, Union

from .ContextInterface import ContextInterface




class Metric(ABC):

    @abstractmethod
    def requires_structure_predictor(self) -> bool:
        raise RuntimeError



    @abstractmethod
    def do(
        self,
        **kwargs,
    ) -> Union[float, Tuple[float]]:
        raise NotImplementedError



    @abstractmethod
    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> Dict[str, float]:
        raise NotImplementedError
