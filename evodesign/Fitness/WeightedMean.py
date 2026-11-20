from __future__ import annotations  # noqa: N999

import numpy as np
import numpy.typing as npt

from ..Evaluation.Evaluator import Evaluator
from .FitnessFunction import FitnessFunction


class WeightedMean(FitnessFunction):
    def __init__(
        self,
        term_calculators: list[Evaluator],
        terms: list[str],
        term_weights: list[float],
        upper_bound: float = 1.0,
    ) -> None:
        super().__init__(upper_bound, terms, term_calculators)
        self.term_weights = term_weights
        assert len(self.term_weights) == len(self.terms)

    def combine(
        self,
        term_values: npt.NDArray[np.float64],
    ) -> float:
        return np.average(term_values, weights=self.term_weights)

    def name(self) -> str:
        return "_".join(
            f"{t.split('.')[-1]}{int(100 * w)}"
            for w, t in zip(self.term_weights, self.terms)
        )
