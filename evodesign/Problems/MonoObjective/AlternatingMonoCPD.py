from .MonoCPD import MonoCPD
from ...Fitness.FitnessFunction import FitnessFunction
from ...Prediction.Predictor import Predictor
from ...Utils.Chain import Chain
from ...Prediction.DirectoryManager import DirectoryManager


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

    def alternate_fitness_fn(self, fitness_fn: FitnessFunction) -> None:
        self.alt_fitness_fn = self.fitness_fn
        self.fitness_fn = fitness_fn
