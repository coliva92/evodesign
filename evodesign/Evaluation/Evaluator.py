from abc import ABC, abstractmethod  # noqa: N999

from .ContextInterface import ContextInterface


class Evaluator(ABC):
    @abstractmethod
    def requires_structure_prediction(self) -> bool:
        raise RuntimeError

    @abstractmethod
    def do(
        self,
        **kwargs,
    ) -> tuple[float, ...]:
        raise NotImplementedError

    @abstractmethod
    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> dict[str, float]:
        raise NotImplementedError
