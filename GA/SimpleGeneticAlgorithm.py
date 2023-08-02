from ..Algorithm import Algorithm
from typing import List, Optional, Dict, Tuple
from ..Fitness import FitnessFunction
from ..Prediction import Predictor
from .Selection import Selection
from .Recombination import Recombination
from .Mutation import Mutation
from .Replacement import GA_Replacement_Generational
from ..Individual import Individual
from evodesign.Statistics import Statistics as Stats
import evodesign.Statistics as Statistics





class SimpleGeneticAlgorithm(Algorithm):

  @classmethod
  def get_name(cls) -> str:
    return 'GA_Simple'



  def __init__(self,
               workspaceName: str,
               targetPdbFilename: str,
               predictor: Predictor,
               fitnessFunction: FitnessFunction,
               populationSize: int,
               numIterations: int,
               selection: Selection,
               recombination: Recombination,
               mutation: Mutation,
               populationFilenames: Optional[List[str]] = None
               ) -> None:
    super().__init__(workspaceName, 
                     targetPdbFilename, 
                     predictor, 
                     fitnessFunction,
                     populationFilenames)
    self._population_size = populationSize
    self._num_iterations = numIterations
    self._selection = selection
    self._recombination = recombination
    self._mutation = mutation
    self._replacement = GA_Replacement_Generational()
    self._statistics = []

  

  def initialize(self) -> List[Individual]:
    m, n = self._sequence_length, self._population_size
    population = [ Individual.new_random(m) for _ in range(n) ]
    self.workspace.save_population(0, population)
    try:
      self._compute_population_fitness(population)
    except RuntimeError as e:
      self.workspace.save_population(0, population)
      raise e
    population = sorted(population)
    self.best_solution = population[-1]
    return population



  def create_next_population(self, 
                             population: List[Individual]
                             ) -> List[Individual]:
    children = self.workspace.restore_children_from_backup()
    if len(children) == 0:
      children = self._evolutionary_step(population)
      self.workspace.backup_children(children)
    try:
      self._compute_population_fitness(children)
    except RuntimeError as e:
      self.workspace.backup_children(children)
      raise e
    self.workspace.delete_children_backup()
    return self._replacement.apply(population, children)



  # TODO: esta funcion, o parte de ella, probablemente se debe mover a la clase 
  # padre
  def run(self, 
          iterationId: int = 0, 
          population: Optional[List[Individual]] = None) -> None:
    if population is None: population = []
    if len(population) == 0:
      iterationId = 0
      population = self.initialize()
      stats = Statistics.compute_statistics(population)
      self._statistics.append(stats)
      self.workspace.save_population(iterationId, population)
      Statistics.save_statistics_to_csv(stats, 
                                        iterationId,
                                        self.workspace.stats_filename)
    if iterationId == 0:
      try:
        is_recovering = self._compute_population_fitness(population)
      except RuntimeError as e:
        self.workspace.save_population(0, population)
        raise e
      if is_recovering:
        population = sorted(population)
        self.best_solution = population[-1]
        stats = Statistics.compute_statistics(population)
        self._statistics.append(stats)
        self.workspace.save_population(iterationId, population)
        Statistics.save_statistics_to_csv(stats, 
                                          iterationId, 
                                          self.workspace.stats_filename)
      else:
        self.best_solution = population[-1]
    else:
      self.best_solution = population[-1]
      self._statistics = Statistics.load_statistics_from_csv(
        self.workspace.stats_filename)
    while True:
      if iterationId == self._num_iterations - 1:
        break
      if population[-1].fitness == self._fitness_fn.upper_bound():
        break
      if self._additional_termination_conditions(iterationId):
        break
      population = self.create_next_population(population)
      iterationId += 1
      stats = Statistics.compute_statistics(population)
      self._statistics.append(stats)
      population, stats = self._update_best_solution(population, stats)
      self.workspace.save_population(iterationId, population)
      Statistics.save_statistics_to_csv(stats, 
                                        iterationId, 
                                        self.workspace.stats_filename)
    Statistics.plot_fitness(self.workspace.stats_filename, 
                            self.workspace.graph_filename,
                            self.workspace.name)



  def _get_params_memento(self) -> dict:
    params = super()._get_params_memento()
    params['fitnessFunction'] = self._fitness_fn.get_name()
    params['populationSize'] = self._population_size
    params['numIterations'] = self._num_iterations
    params['selection'] = self._selection.get_name()
    params['selectionParams'] = self._selection.get_params_memento()
    params['recombination'] = self._recombination.get_name()
    params['recombinationParams'] = self._recombination.get_params_memento()
    params['mutation'] = self._mutation.get_name()
    params['mutationParams'] = self._mutation.get_params_memento()
    return params
  


  def _evolutionary_step(self, 
                         population: List[Individual]
                         ) -> List[Individual]:
    next_population = []
    while len(next_population) < len(population):
      parents = self._selection.apply(population)
      children = self._recombination.apply(parents)
      self._mutation.apply(children)
      next_population += children
    return next_population
  

  
  def _compute_population_fitness(self, 
                                  population: List[Individual]
                                  ) -> bool:
    was_some_fitness_computed = False
    for individual in population:
      if individual.fitness is None:
        was_some_fitness_computed = True
        filename = individual.get_pdb_filepath(self.workspace.pdbs_folder)
        individual.update_fitness(self._fitness_fn, 
                                  self._predictor, 
                                  self._reference_backbone,
                                  filename)
    return was_some_fitness_computed



  def _update_best_solution(self, 
                            population: List[Individual], 
                            stats: Stats
                            ) -> Tuple[List[Individual], Stats]:
    if population[-1] > self.best_solution: 
      self.best_solution = population[-1]
    if population[-1] < self.best_solution:
      population = population[1:] + [ self.best_solution ]
    stats.best_sequence_fitness = self.best_solution.fitness
    stats.best_sequence = f'"{self.best_solution.sequence}"'
    return population, stats



  def _additional_termination_conditions(self, iterationId: int) -> bool:
    return False
