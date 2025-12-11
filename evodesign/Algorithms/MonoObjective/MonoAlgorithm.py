from ..Algorithm import Algorithm
from ...Problems.MonoObjective.MonoCPD import MonoCPD
from ...Utils.Chain import Chain
from ...Prediction.Predictor import Predictor
from ...Prediction.DirectoryManager import DirectoryManager
from ...Fitness.FitnessFunction import FitnessFunction


class MonoAlgorithm(Algorithm):

    def __init__(
        self,
        max_generations: int,
        population_size: int,
        fitness_fn: FitnessFunction,
        predictor: Predictor,
    ):
        super().__init__(max_generations, population_size, predictor)
        self.fitness_fn = fitness_fn

    def _create_problem(
        self,
        ref_chain: Chain,
        predictor_directory: DirectoryManager,
    ) -> MonoCPD:
        return MonoCPD(
            ref_chain,
            self.fitness_fn,
            self.predictor,
            predictor_directory,
        )

    def num_terms(self) -> int:
        return self.fitness_fn.num_terms()
