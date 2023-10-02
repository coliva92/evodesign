from ..Algorithm import Algorithm
from typing import Optional, Callable, List
from ..Fitness import FitnessFunction, Fitness_RmsdGdtEnergyScore
from ..Prediction import Predictor
from .Selection import Selection
from .Recombination import Recombination
from .Mutation import Mutation
from .Replacement import GA_Replacement_Generational
from ..Population import Population
from ..Statistics import Statistics
from .Termination import DiversityLowerBoundReached
from ..Individual import Individual





class SimpleGeneticAlgorithm(Algorithm):

  @classmethod
  def name(cls) -> str:
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
               mutation: Mutation,
               elitismSize: int
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
    self._elitismSize = elitismSize
    self._replacement = GA_Replacement_Generational()
    self._terminators = [ DiversityLowerBoundReached() ]
    if isinstance(self._fitness_fn, Fitness_RmsdGdtEnergyScore):
      self._fitness_fn._metric_calculators['energyScore'].set_pdbs_folder(
        self.workspace.pdbs_folder)
  


  def _params_json(self) -> dict:
    a = super()._params_json()
    b = {
      'populationSize': self._population_size,
      'numIterations': self._num_iterations,
      'elitismSize': self._elitismSize,
      'selection': self._selection.name(),
      'selectionParams': self._selection.params_json(),
      'recombination': self._recombination.name(),
      'recombinationParams': self._recombination.params_json(),
      'mutation': self._mutation.name(),
      'mutationParams': self._mutation.params_json()
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
    stats = None
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
      print(f'{population.iteration_id:04d} / {self._num_iterations:04d} ' + \
            f'{self.best_solution.fitness:.5f} ' + \
            f'{stats.diversity:0.5f}',
            flush=True)
    while True:
      if population.iteration_id == self._num_iterations:
        break
      if self.best_solution.fitness >= self._fitness_fn.upper_bound():
        break
      if sum([ term(population, stats) for term in self._terminators ]):
        break
      next_population = self.next_population(population)
      self.workspace.save_population(population)
      stats = Statistics.new_from_population(next_population)
      self._update_best_solution(next_population, population)
      population = next_population
      self.workspace.save_statistics(stats, self.best_solution)
      print(f'{population.iteration_id:04d} / {self._num_iterations:04d} ' + \
            f'{self.best_solution.fitness:.5f} ' + \
            f'{stats.diversity:0.5f}',
            flush=True)



  def _update_best_solution(self, 
                            next_population: Population,
                            old_population: Population
                            ) -> None:
    top = self._merge(next_population[-self._elitismSize:],
                      old_population[-self._elitismSize:])
    if top[-1] > self.best_solution: 
      self.best_solution = top[-1]
  


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
  


  def _merge(self, 
             a: List[Individual], 
             b: List[Individual]
             ) -> List[Individual]:
    merged, i, j = [], 0, 0
    while True:
      if i == len(a):
        merged += b[j:]
        break
      if j == len(b):
        merged += a[i:]
        break
      if a[i] <= b[j]:
        merged.append(a[i])
        i += 1
      else:
        merged.append(b[j])
        j += 1
    return merged
    