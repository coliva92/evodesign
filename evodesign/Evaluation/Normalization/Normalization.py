from abc import ABC, abstractmethod  # noqa: N999


class Normalization(ABC):
    @abstractmethod
    def do(self, x: float, **kwargs) -> float:
        raise NotImplementedError
