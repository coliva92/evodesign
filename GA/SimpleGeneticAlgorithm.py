from ..Algorithm import Algorithm
from typing import List, Optional, Dict, Tuple
from ..Fitness import FitnessFunction
from ..Prediction import Predictor
from .Selection import Selection
from .Recombination import Recombination
from .Mutation import Mutation
from ..Population import Individual
import evodesign.Statistics as Statistics





class SimpleGeneticAlgorithm(Algorithm):

  def __init__(self,
               workspaceName: str,
               targetPdbFilename: str,
               predictor: Predictor,
               fitnessFunction: FitnessFunction,
               populationSize: int,
               numIterations: int,
               selection: Selection,
               recombination: Recombination,
               mutation: Mutation
               ) -> None:
    super().__init__(workspaceName, targetPdbFilename, predictor)
    self._population_size = populationSize
    self._num_iterations = numIterations
    self._fitness_fn = fitnessFunction
    self._selection = selection
    self._recombination = recombination
    self._mutation = mutation
    self._best_solution = None



  @classmethod
  def get_type(cls) -> str:
    return 'GA_Simple'



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

  

  def initialize(self) -> List[Individual]:
    m, n = self._sequence_length, self._population_size
    population = [ Individual.random(m) for _ in range(n) ]
    self.workspace.save_population(0, population)
    self._compute_initial_population_fitness(population)
    return sorted(population)



  def step(self, population: List[Individual]) -> List[Individual]:
    if self.workspace.has_backed_up_children():
      new_population = self.workspace.restore_children_from_backup()
    else:
      new_population = []
      while len(new_population) < len(population):
        parents = self._selection.apply(population)
        children = self._recombination.apply(parents)
        self._mutation.apply(children)
        new_population += children
      self.workspace.backup_children(new_population)
    for i, child in enumerate(new_population):
      if child.fitness == None:
        filename = self._get_pdb_filename(child)
        self._fitness_fn.apply(child, 
                              self._predictor, 
                              self._reference_backbone, 
                              filename)
        self.workspace.update_children_backup(i, child)
    self.workspace.delete_children_backup()
    return sorted(new_population)



  # TODO: esta funcion, o parte de ella, probablemente se debe mover a la clase 
  # padre
  def run(self, 
          iterationId: int = 0, 
          population: Optional[List[Individual]] = None) -> None:
    if population == None: population = []
    is_recovering = False
    if len(population) == 0:
      iterationId = 0
      population = self.initialize()
      stats = Statistics.compute_statistics(population)
      self._best_solution = population[-1]
      self.workspace.save_population(iterationId, population)
      Statistics.save_statistics_to_csv_file(stats, 
                                             iterationId, 
                                             self.workspace.stats_filename)
    else:
      self._recovering_statistics()
      is_recovering = self._compute_initial_population_fitness(population)
    if is_recovering:
      population = sorted(population)
      stats = Statistics.compute_statistics(population)
      population, stats = self._update_best_solution(population, stats)
      self.workspace.save_population(iterationId, population)
      Statistics.save_statistics_to_csv_file(stats, 
                                             iterationId, 
                                             self.workspace.stats_filename)
    while True:
      if iterationId == self._num_iterations - 1:
        break
      if population[-1].fitness == self._fitness_fn.upper_bound():
        break
      population = self.step(population)
      iterationId += 1
      stats = Statistics.compute_statistics(population)
      population, stats = self._update_best_solution(population, stats)
      self.workspace.save_population(iterationId, population)
      Statistics.save_statistics_to_csv_file(stats, 
                                             iterationId, 
                                             self.workspace.stats_filename)
    Statistics.plot_fitness_over_iterations(self.workspace.stats_filename, 
                                            self.workspace.graph_filename,
                                            self.workspace.name)



  def _compute_initial_population_fitness(self, 
                                          population: List[Individual]
                                          ) -> bool:
    computed_some_fitness = False
    for individual in population:
      if individual.fitness == None:
        computed_some_fitness = True
        try:
          filename = self._get_pdb_filename(individual)
          self._fitness_fn.apply(individual, 
                                 self._predictor, 
                                 self._reference_backbone,
                                 filename)
        except RuntimeError as e:
          self.workspace.save_population(0, population)
          raise e
    return computed_some_fitness



  def _recovering_statistics(self) -> None:
    temp = Statistics.load_statistics_from_csv_file(
      self.workspace.stats_filename)
    if len(temp) > 0:
      self._best_solution = Individual(temp[-1]['best_sequence'],
                                      temp[-1]['best_sequence_id'],
                                      temp[-1]['best_sequence_fitness'])



  def _update_best_solution(self, 
                            population: List[Individual], 
                            stats: Dict[str, float]
                            ) -> Tuple[List[Individual], Dict[str, float]]:
    if self._best_solution == None or population[-1] > self._best_solution: 
      self._best_solution = population[-1]
    elif population[-1] != self._best_solution:
      population = population[1:] + [ self._best_solution ]
    stats['best_sequence_id'] = self._best_solution.id
    stats['best_sequence_fitness'] = self._best_solution.fitness
    stats['best_sequence'] = f'"{self._best_solution.sequence}"'
    return population, stats
  