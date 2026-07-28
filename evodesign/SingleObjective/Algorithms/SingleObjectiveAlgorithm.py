from typing import Optional

import numpy as np
import numpy.typing as npt

from ...AlgorithmFactory import AlgorithmFactory
from ...Chemistry.Chain import Chain
from ...Fitness.FitnessFunction import FitnessFunction
from ...Prediction.DirectoryManager import DirectoryManager
from ...SingleObjective.Problems.SingleObjectiveCPD import SingleObjectiveCPD


class SingleObjectiveAlgorithm(AlgorithmFactory):

    def __init__(
        self,
        fitness_fn: FitnessFunction,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.fitness_fn = fitness_fn
        return

    def create_problem(
        self,
        ref_chain: Chain,
        predictor_directory: DirectoryManager,
        aa_profile: Optional[npt.NDArray[np.float64]] = None,
    ) -> SingleObjectiveCPD:
        return SingleObjectiveCPD(
            ref_chain=ref_chain,
            fitness_fn=self.fitness_fn,
            predictor=self.predictor,
            predictor_directory=predictor_directory,
            aa_profile=aa_profile,
        )

    def num_terms(self) -> int:
        return self.fitness_fn.num_terms()
