from typing import Optional

from pymoo.algorithms.soo.nonconvex.ga import GA

from ....GA.Crossover.Crossover import Crossover
from ....GA.Crossover.UniformCrossover import UniformCrossover
from ....GA.Mutation.Mutation import Mutation
from ....GA.Mutation.RandomResetting import RandomResetting
from ....GA.Replacement.PyMOO.Generational import \
    Generational as GenerationalReplacement
from ....GA.Selection.Selection import Selection
from ....GA.Selection.Tournament import Tournament
from ....GA.Termination.MaximumDiversityLoss import MaximumDiversityLoss
from ..SingleObjectiveAlgorithm import SingleObjectiveAlgorithm




class Generational(SingleObjectiveAlgorithm):

    def __init__(
        self,
        diversity_loss_tol: float = 0.95,
        sample_size: Optional[int] = 30,
        selection: Selection = Tournament(),
        crossover: Crossover = UniformCrossover(),
        mutation: Mutation = RandomResetting(),
        **kwargs
    ) -> None:
        super().__init__(**kwargs)
        self.selection = selection
        self.crossover = crossover
        self.mutation = mutation
        self.diversity_loss_tol = diversity_loss_tol
        self.sample_size = sample_size
        return



    def create_algorithm(self) -> GA:
        replacement = GenerationalReplacement()
        termination = MaximumDiversityLoss(
            self.max_generations, self.diversity_loss_tol, self.sample_size
        )
        algorithm = GA(
            pop_size=self.population_size,
            n_offsprings=self.population_size,
            sampling=self._sampling,
            selection=self.selection._pymoo_selection,
            crossover=self.crossover._pymoo_crossover,
            mutation=self.mutation._pymoo_mutation,
            survival=replacement,
            eliminate_duplicates=False,
        )
        algorithm.termination = termination
        return algorithm
