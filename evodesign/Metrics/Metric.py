from abc import ABC, abstractmethod
from ..RetrievableSettings import RetrievableSettings
from .ContextInterface import ContextInterface
from typing import Union, Tuple, Dict


class Metric(RetrievableSettings, ABC):

    @abstractmethod
    def uses_predictor(self) -> bool:
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
