from .MonoCPD import MonoCPD
from ...Fitness.FitnessFunction import FitnessFunction
from ...Prediction.Predictor import Predictor
from ...Utils.Chain import Chain
from ...Prediction.DirectoryManager import DirectoryManager
from typing import List
import numpy as np
import numpy.typing as npt


class AlternatingMonoCPD(MonoCPD):

    def __init__(
        self,
        ref_chain: Chain,
        fitness_fn: FitnessFunction,
        alt_fitness_fn: FitnessFunction,
        predictor: Predictor,
        predictor_directory: DirectoryManager,
    ):
        super().__init__(ref_chain, fitness_fn, predictor, predictor_directory)
        self.alt_fitness_fn = alt_fitness_fn
        self._curr_fn_idx = 0
        self._functions_list = [self.fitness_fn, self.alt_fitness_fn]
        return

    def _compute_term_values(
        self,
        model_chains: List[Chain],
    ) -> npt.NDArray[np.float64]:
        f = self._functions_list[self._curr_fn_idx]
        g = self._functions_list[int(not self._curr_fn_idx)]
        term_values = []
        for i in range(len(self._functions_list)):
            if i == self._curr_fn_idx:
                term_values.extend(
                    [f.do(model_chain, self.ref_chain) for model_chain in model_chains]
                )
            else:
                term_values.extend([None for _ in range(g.num_terms())])
        return np.array(term_values)

    def alternate_fitness_fn(self) -> None:
        self._curr_fn_idx = int(not self._curr_fn_idx)
        return
