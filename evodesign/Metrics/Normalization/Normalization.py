from abc import ABC, abstractmethod
from ...RetrievableSettings import RetrievableSettings


class Normalization(RetrievableSettings, ABC):

    @abstractmethod
    def do(self, x: float) -> float:
        raise NotImplementedError
