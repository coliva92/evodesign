from pymoo.algorithms.soo.nonconvex.ga import GA
from pymoo.termination.max_eval import MaximumFunctionCallTermination

from ....Callbacks.CallbackCollection import CallbackCollection
from ....Callbacks.PopulationRestarter import PopulationRestarter
from ....Callbacks.StorageManager import StorageManager
from ....GA.Crossover.Crossover import Crossover
from ....GA.Crossover.UniformCrossover import UniformCrossover
from ....GA.Mutation.RandomResetting import RandomResetting
from ....GA.Replacement.PyMOO.GenerationalElitism import GenerationalElitism
from ....GA.Selection.Selection import Selection
from ....GA.Selection.Tournament import Tournament
from ..SingleObjectiveAlgorithm import SingleObjectiveAlgorithm




class Micro(SingleObjectiveAlgorithm):

    def __init__(
        self,
        max_fitness_fn_evals: int,
        diversity_loss_tol: float = 0.95,
        selection: Selection = Tournament(),
        crossover: Crossover = UniformCrossover(),
        **kwargs
    ) -> None:
        max_generations = (max_fitness_fn_evals // kwargs["popultion_size"]) + 1
        super().__init__(max_generations=max_generations, **kwargs)
        self.max_fitness_fn_evals = max_fitness_fn_evals
        self.diversity_loss_tol = diversity_loss_tol
        self.selection = selection
        self.crossover = crossover
        return



    def create_algorithm(self) -> GA:
        mutation = RandomResetting(sequence_mutation_prob=0, residue_mutation_prob=0)
        replacement = GenerationalElitism()
        termination = MaximumFunctionCallTermination(self.max_fitness_fn_evals)
        algorithm = GA(
            pop_size=self.population_size,
            n_offsprings=self.population_size,
            sampling=self._sampling,
            selection=self.selection._pymoo_selection,
            crossover=self.crossover._pymoo_crossover,
            mutation=mutation._pymoo_mutation,
            survival=replacement,
            eliminate_duplicates=False,
        )
        algorithm.termination = termination
        return algorithm



    def create_callbacks(
        self,
        storage: StorageManager,
    ) -> CallbackCollection:
        callbacks = super().create_callbacks(storage)
        callbacks.append(PopulationRestarter(self.diversity_loss_tol))
        return callbacks
