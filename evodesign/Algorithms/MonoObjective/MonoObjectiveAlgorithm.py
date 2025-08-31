from ..Algorithm import Algorithm
from .MonoObjectiveCPD import MonoObjectiveCPD
from ...Utils.Chain import Chain
from ...Prediction.Predictor import Predictor
from ...Fitness.FitnessFunction import FitnessFunction


class MonoObjectiveAlgorithm(Algorithm):

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
        prediction_pdb_path: str,
    ) -> MonoObjectiveCPD:
        return MonoObjectiveCPD(
            self.predictor,
            self.fitness_fn,
            ref_chain,
            prediction_pdb_path,
        )
    
    def num_terms(self) -> int:
        return self.fitness_fn.num_terms()
