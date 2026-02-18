from ..Algorithm import Algorithm
from ...Problems.MonoObjective.MonoCPD import MonoCPD
from ...Chemistry.Chain import Chain
from ...Prediction.Predictor import Predictor
from ...Prediction.DirectoryManager import DirectoryManager
from ...Fitness.FitnessFunction import FitnessFunction
from typing import Optional
import numpy as np
import numpy.typing as npt


class MonoAlgorithm(Algorithm):

    def __init__(
        self,
        max_generations: int,
        population_size: int,
        predictor: Predictor,
        fitness_fn: FitnessFunction,
    ):
        super().__init__(max_generations, population_size, predictor)
        self.fitness_fn = fitness_fn

    def _create_problem(
        self,
        ref_chain: Chain,
        predictor_directory: DirectoryManager,
        aa_profile: Optional[npt.NDArray[np.float64]] = None,
    ) -> MonoCPD:
        return MonoCPD(
            ref_chain,
            self.fitness_fn,
            self.predictor,
            predictor_directory,
            aa_profile
        )

    def num_terms(self) -> int:
        return self.fitness_fn.num_terms()
