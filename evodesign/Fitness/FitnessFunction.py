from __future__ import annotations  # noqa: N999

from abc import ABC, abstractmethod

import numpy as np
import numpy.typing as npt

from ..Evaluation.Context import Context
from ..Evaluation.Evaluator import Evaluator
from ..Protein.Chain import Chain
from ..Settings import get_namespace


class FitnessFunction(ABC):
    def __init__(
        self,
        upper_bound: float,
        terms: list[str],
        evaluators: list[Evaluator],
    ) -> None:
        super().__init__()
        self.upper_bound = upper_bound
        self.terms = terms
        self.evaluators = evaluators
        self._evaluators = {get_namespace(calc): calc for calc in self.evaluators}

    def do(
        self,
        model_chain: Chain,
        ref_chain: Chain,
        persistent_params: dict,
        **kwargs,
    ) -> npt.NDArray[np.float64]:
        context = Context(model_chain, ref_chain, self._evaluators, persistent_params, **kwargs)
        term_values = np.array(
            [context.get_component_value(name) for name in self.terms]
        )
        return np.concatenate((np.array([self.combine(term_values)]), term_values))

    @abstractmethod
    def combine(
        self,
        term_values: npt.NDArray[np.float64],
    ) -> float:
        raise NotImplementedError

    @abstractmethod
    def name(self) -> str:
        raise NotImplementedError

    def num_terms(self) -> int:
        return len(self.terms)

    def requires_structure_prediction(self) -> bool:
        active_terms = {".".join(term.split(".")[:-1]) for term in self.terms}
        return any(
            self._evaluators[name].requires_structure_prediction()
            for name in active_terms
        )
