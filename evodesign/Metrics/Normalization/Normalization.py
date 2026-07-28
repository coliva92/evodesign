from abc import ABC, abstractmethod




class Normalization(ABC):

    @abstractmethod
    def do(self, x: float) -> float:
        raise NotImplementedError
