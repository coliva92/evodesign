from ..MonoAlgorithm import MonoAlgorithm
from ....Chemistry.Chain import Chain
from ....Fitness.FitnessFunction import FitnessFunction
from ....Prediction.Predictor import Predictor
from ....GA.Selection.Selection import Selection
from ....GA.Selection.Tournament import Tournament
from ....GA.Crossover.Crossover import Crossover
from ....GA.Crossover.UniformCrossover import UniformCrossover
from ....GA.Mutation.RandomResetting import RandomResetting
from ....GA.Replacement.PyMOO.GenerationalElitism import GenerationalElitism
from ....Callbacks.StorageManager import StorageManager
from ....Callbacks.PopulationRestarter import PopulationRestarter
from pymoo.algorithms.soo.nonconvex.ga import GA
from pymoo.core.callback import Callback
import numpy as np
import numpy.typing as npt
from typing import Optional


class Micro(MonoAlgorithm):

    def __init__(
        self,
        max_fitness_fn_evals: int,
        population_size: int,
        predictor: Predictor,
        fitness_fn: FitnessFunction,
        diversity_loss_tol: float = 0.95,
        selection: Selection = Tournament(),
        crossover: Crossover = UniformCrossover(),
    ):
        max_generations = np.ceil(max_fitness_fn_evals / population_size, dtype=np.int64)
        super().__init__(max_generations, population_size, predictor, fitness_fn)
        self.max_fitness_fn_evals = max_fitness_fn_evals
        self.diversity_loss_tol = diversity_loss_tol
        self.selection = selection
        self.crossover = crossover
        self._replacement = GenerationalElitism()

    def _create_algorithm(self) -> GA:
        return GA(
            pop_size=self.population_size,
            n_offsprings=self.population_size - 1,
            sampling=self._sampling,
            selection=self.selection._pymoo_selection,
            crossover=self.crossover._pymoo_crossover,
            mutation=RandomResetting(
                sequence_mutation_prob=0, residue_mutation_prob=0
            )._pymoo_mutation,
            survival=self._replacement,
            eliminate_duplicates=False,
            termination=("n_eval", self.max_fitness_fn_evals),
        )

    def _callbacks_chain(
        self,
        storage: StorageManager,
    ) -> Callback:
        return PopulationRestarter(self.diversity_loss_tol, storage)

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
            )
        return super().run(ref_chain, storage, **kwargs)
