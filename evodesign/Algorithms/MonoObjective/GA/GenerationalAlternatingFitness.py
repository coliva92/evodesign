from .Generational import Generational
from pymoo.core.callback import Callback
from ....Utils.Chain import Chain
from ....Prediction.DirectoryManager import DirectoryManager
from ....Problems.MonoObjective.AlternatingMonoCPD import AlternatingMonoCPD
from ....Callbacks.StorageManager import StorageManager
from ....Callbacks.MonoFitnessFnAlternator import MonoFitnessFnAlternator
from ....Fitness.FitnessFunction import FitnessFunction
from ....Prediction.Predictor import Predictor
from ....GA.Selection.Selection import Selection
from ....GA.Selection.Tournament import Tournament
from ....GA.Crossover.Crossover import Crossover
from ....GA.Crossover.UniformCrossover import UniformCrossover
from ....GA.Mutation.Mutation import Mutation
from ....GA.Mutation.RandomResetting import RandomResetting
from typing import Optional
import numpy as np
import numpy.typing as npt


class GenerationalAlternatingFitness(Generational):

    def __init__(
        self,
        max_generations: int,
        population_size: int,
        predictor: Predictor,
        fitness_fn: FitnessFunction,
        alt_fitness_fn: FitnessFunction,
        selection: Selection = Tournament(),
        crossover: Crossover = UniformCrossover(),
        mutation: Mutation = RandomResetting(),
        alt_fitness_fn_every_nth_generation: int = 10,
    ):
        super().__init__(
            max_generations,
            population_size,
            predictor,
            fitness_fn,
            selection,
            crossover,
            mutation,
        )
        self.alt_fitness_fn = alt_fitness_fn
        self.alt_fitness_fn_every_nth_generation = alt_fitness_fn_every_nth_generation
        return

    def num_terms(self):
        return self.fitness_fn.num_terms() + self.alt_fitness_fn.num_terms()

    def _create_problem(
        self,
        ref_chain: Chain,
        predictor_directory: DirectoryManager,
        aa_profile: Optional[npt.NDArray[np.float64]] = None,
    ) -> AlternatingMonoCPD:
        return AlternatingMonoCPD(
            ref_chain,
            self.fitness_fn,
            self.alt_fitness_fn,
            self.predictor,
            predictor_directory,
            aa_profile,
        )

    def _callbacks_chain(
        self,
        storage: StorageManager,
    ) -> Callback:
        return MonoFitnessFnAlternator(
            self.alt_fitness_fn_every_nth_generation,
            storage,
        )
