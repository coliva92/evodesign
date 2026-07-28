from typing import List

import numpy as np
import numpy.typing as npt

from ...Chemistry.Chain import Chain
from ...Fitness.FitnessFunction import FitnessFunction
from .SingleObjectiveCPD import SingleObjectiveCPD


class AlternatingCPD(SingleObjectiveCPD):

    _NUM_FITNESS_FUNCS = 2

    def __init__(self, alt_fitness_fn: FitnessFunction, **kwargs):
        super().__init__(**kwargs)
        self.alt_fitness_fn = alt_fitness_fn
        self.alt_archive = {}
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
            row = [term_values[0]]
            for i in range(self._NUM_FITNESS_FUNCS):
                if i == self._curr_fn_idx:
                    row.extend(term_values[1:].tolist())
                else:
                    row.extend([0 for _ in range(g.num_terms())])
            terms_matrix.append(row)
        return np.array(terms_matrix)

    def alternate_fitness_fn(self) -> None:
        self._curr_fn_idx = int(not self._curr_fn_idx)
        tmp_f = self.fitness_fn
        tmp_arch = self.archive
        self.fitness_fn = self.alt_fitness_fn
        self.archive = self.alt_archive
        self.alt_fitness_fn = tmp_f
        self.alt_archive = tmp_arch
        return
