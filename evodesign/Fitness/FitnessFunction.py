from abc import ABC, abstractmethod
from ..RetrievableSettings import RetrievableSettings
from ..Metrics.Metric import Metric
from ..Metrics.Context import Context
from ..Utils.Chain import Chain
import numpy as np
import numpy.typing as npt
from typing import List


class FitnessFunction(RetrievableSettings, ABC):

    def __init__(
        self,
        upper_bound: float,
        terms: List[str],
        term_calculators: List[Metric],
    ) -> None:
        super().__init__()
        self.upper_bound = upper_bound
        self.terms = terms
        self.term_calculators = term_calculators
        self._term_calculators = {
            calc._class_name(): calc for calc in self.term_calculators
        }

    def do(
        self,
        model_chain: Chain,
        ref_chain: Chain,
        **kwargs,
    ) -> npt.NDArray[np.float64]:
        context = Context(model_chain, ref_chain, self._term_calculators, **kwargs)
        term_values = np.array(
            [context.get_component_value(name) for name in self.terms]
        )
        component_values = [
            component_value
            for term_name in self.terms
            for _, component_value in context.get_metric_components(
                ".".join(term_name.split(".")[:-1])
            ).items()
        ]
        return np.array([self.combine(term_values), *component_values])

    @abstractmethod
    def combine(
        self,
        term_values: npt.NDArray[np.float64],
    ) -> float:
        raise NotImplementedError

    def num_terms(self) -> int:
        return len(self.terms)
