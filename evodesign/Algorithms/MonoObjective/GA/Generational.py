from ..MonoObjectiveAlgorithm import MonoObjectiveAlgorithm
from ....Utils.Chain import Chain
from ....Fitness.FitnessFunction import FitnessFunction
from ....Prediction.Predictor import Predictor
from ....GA.Selection.Selection import Selection
from ....GA.Selection.Tournament import Tournament
from ....GA.Crossover.Crossover import Crossover
from ....GA.Crossover.UniformCrossover import UniformCrossover
from ....GA.Mutation.Mutation import Mutation
from ....GA.Mutation.RandomResetting import RandomResetting
from pymoo.termination.max_gen import MaximumGenerationTermination
from ....GA.Replacement.PyMOO.Generational import (
    Generational as GenerationalReplacement,
)
from ....Utils.StorageManager import StorageManager
from pymoo.algorithms.soo.nonconvex.ga import GA
from typing import Optional


class Generational(MonoObjectiveAlgorithm):

    def __init__(
        self,
        max_generations: int,
        population_size: int,
        predictor: Predictor,
        fitness_fn: FitnessFunction,
        selection: Selection = Tournament(),
        crossover: Crossover = UniformCrossover(),
        mutation: Mutation = RandomResetting(),
    ):
        super().__init__(max_generations, population_size, predictor, fitness_fn)
        self.selection = selection
        self.crossover = crossover
        self.mutation = mutation
        self._termination = MaximumGenerationTermination(self.max_generations)
        self._replacement = GenerationalReplacement()

    def _create_algorithm(self) -> GA:
        return GA(
            pop_size=self.population_size,
            n_offsprings=self.population_size,
            sampling=self._sampling,
            selection=self.selection._pymoo_selection,
            crossover=self.crossover._pymoo_crossover,
            mutation=self.mutation._pymoo_mutation,
            survival=self._replacement,
            eliminate_duplicates=False,
        )

    def run(
        self,
        ref_chain: Chain,
        storage: StorageManager,
        **kwargs,
    ):
        if self._algorithm is None or self._algorithm.termination is None:
            return super().run(
                ref_chain, storage, **{**{"termination": self._termination}, **kwargs}
            )
        return super().run(ref_chain, storage, **kwargs)
