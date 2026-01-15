from .MonoCPD import MonoCPD
from ...Fitness.FitnessFunction import FitnessFunction
from ...Prediction.Predictor import Predictor
from ...Utils.Chain import Chain
from ...Prediction.DirectoryManager import DirectoryManager
from typing import List
import numpy as np
import numpy.typing as npt


class AlternatingMonoCPD(MonoCPD):

    _NUM_FITNESS_FUNCS = 2

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
        return

    def _compute_term_values(
        self,
        model_chains: List[Chain],
    ) -> npt.NDArray[np.float64]:
        f = self.fitness_fn
        g = self.alt_fitness_fn
        terms_matrix = []
        for model_chain in model_chains:
            term_values = f.do(model_chain, self.ref_chain)
            row = [ term_values[0] ]
            for i in range(self._NUM_FITNESS_FUNCS):
                if i == self._curr_fn_idx:
                    row.extend(term_values[1:].tolist())
                else:
                    row.extend([ 0 for _ in range(g.num_terms()) ])
            terms_matrix.append(row)            
        return np.array(terms_matrix)

    def alternate_fitness_fn(self) -> None:
        self._curr_fn_idx = int(not self._curr_fn_idx)
        tmp = self.fitness_fn
        self.fitness_fn = self.alt_fitness_fn
        self.alt_fitness_fn = tmp
        return
