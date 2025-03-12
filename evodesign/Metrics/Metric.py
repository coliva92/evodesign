from abc import ABC, abstractmethod
from ..RetrievableSettings import RetrievableSettings


class Metric(RetrievableSettings, ABC):

    @abstractmethod
    def do(self, **kwargs) -> float:
        raise NotImplementedError
