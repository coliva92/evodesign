from ..MonoAlgorithm import MonoAlgorithm
from ....Utils.Chain import Chain
from ....Fitness.FitnessFunction import FitnessFunction
from ....Prediction.Predictor import Predictor
from ....GA.Selection.Selection import Selection
from ....GA.Selection.Tournament import Tournament
from ....GA.Crossover.Crossover import Crossover
from ....GA.Crossover.UniformCrossover import UniformCrossover
from ....GA.Mutation.Mutation import Mutation
from ....GA.Mutation.RandomResetting import RandomResetting
from ....GA.Termination.MaximumDiversityLoss import MaximumDiversityLoss
from ....GA.Replacement.PyMOO.Generational import (
    Generational as GenerationalReplacement,
)
from ....Callbacks.StorageManager import StorageManager
from pymoo.algorithms.soo.nonconvex.ga import GA
import numpy as np
import numpy.typing as npt
from typing import Optional


class Generational(MonoAlgorithm):

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
        self._termination = MaximumDiversityLoss(
            self.max_generations, max_similarity=0.95, sample_size=30
        )
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
        aa_profile: Optional[npt.NDArray[np.float64]] = None,
        **kwargs,
    ):
        if self._algorithm is None or self._algorithm.termination is None:
            return super().run(
                ref_chain,
                storage,
                aa_profile,
                **{**{"termination": self._termination}, **kwargs},
            )
        return super().run(ref_chain, storage, **kwargs)
