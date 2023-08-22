from ..Algorithm import Algorithm
from typing import Optional, Callable
from ..Fitness import FitnessFunction
from ..Prediction import Predictor
from .Selection import Selection
from .Recombination import Recombination
from .Mutation import Mutation
from .Replacement import GA_Replacement_Generational
from ..Population import Population
from ..Statistics import Statistics





class SimpleGeneticAlgorithm(Algorithm):

  @classmethod
  def get_name(cls) -> str:
    return 'GA_Simple'



  def __init__(self,
               workspaceRoot: str,
               targetPdbFilename: str,
               predictor: Predictor,
               fitnessFunction: FitnessFunction,
               populationSize: int,
               numIterations: int,
               selection: Selection,
               recombination: Recombination,
               mutation: Mutation
               ) -> None:
    super().__init__(workspaceRoot, 
                     targetPdbFilename, 
                     predictor, 
                     fitnessFunction)
    self._population_size = populationSize
    self._num_iterations = numIterations
    self._selection = selection
    self._recombination = recombination
    self._mutation = mutation
    self._replacement = GA_Replacement_Generational()
    self._terminators = []
  


  def _params_json(self) -> dict:
    a = super()._params_json()
    b = {
      'populationSize': self._population_size,
      'numIterations': self._num_iterations,
      'selection': self._selection.get_name(),
      'selectionParams': self._selection.as_json(),
      'recombination': self._recombination.get_name(),
      'recombinationParams': self._recombination.as_json(),
      'mutation': self._mutation.get_name(),
      'mutationParams': self._mutation.as_json()
    }
    return { **a, **b }
  


  def _evolutionary_step(self, population: Population) -> Population:
    next_population = Population()
    while len(next_population) < len(population):
      parents = self._selection(population)
      children = self._recombination(parents)
      self._mutation(children)
      next_population += children
    return next_population



  def next_population(self, population: Population) -> Population:
    children = self.workspace.restore_children_from_backup()
    if not children:
      children = self._evolutionary_step(population)
      self.workspace.backup_children(children)
      self.workspace.save_rng_settings()
    self._update_fitness(children, self.workspace.backup_children)
    self.workspace.delete_children_backup()
    self.workspace.save_rng_settings()
    return self._replacement(population, children)



  def __call__(self, population: Optional[Population] = None) -> None:
    if population is None: population = Population()
    if not population:
      population = Population.new_random(self._population_size, 
                                         self._sequence_length)
      self.workspace.save_population(population)
      self.workspace.save_rng_settings()
    is_recovering = self._update_fitness(population, 
                                         self.workspace.save_population)
    self.best_solution = population[-1]
    if is_recovering:
      stats = Statistics.new_from_population(population)
      self.workspace.save_population(population)
      self.workspace.save_statistics(stats, self.best_solution)
      self.workspace.save_rng_settings()
    while True:
      if population.iteration_id == self._num_iterations:
        break
      if self.best_solution.fitness >= self._fitness_fn.upper_bound():
        break
      if sum([ term(population, stats) for term in self._terminators ]):
        break
      population = self.next_population(population)
      stats = Statistics.new_from_population(population)
      self._update_best_solution(population)
      self.workspace.save_population(population)
      self.workspace.save_statistics(stats, self.best_solution)



  def _update_best_solution(self, population: Population) -> None:
    if population[-1] > self.best_solution: 
      self.best_solution = population[-1]
    elif population[-1].sequence != self.best_solution.sequence:
      population.individuals = population[1:] + [ self.best_solution ]
  


  def _update_fitness(self, 
                      population: Population, 
                      saveFunction: Callable[[Population], None]
                      ) -> bool:
    try:
      result = population.update_fitness(self._fitness_fn, 
                                         self._predictor, 
                                         self._reference_backbone, 
                                         self.workspace.pdbs_folder)
    except BaseException as e:
      saveFunction(population)
      self.workspace.save_rng_settings()
      raise e
    return result
    