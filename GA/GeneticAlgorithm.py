from ..Algorithm import Algorithm
from typing import List, Optional
from ..Fitness import FitnessFunction
from ..Prediction import Predictor
from .Selection import Selection
from .Recombination import Recombination
from .Mutation import Mutation
from .Replacement import Replacement
from ..Population import Individual
import evodesign.Statistics as Statistics





class GeneticAlgorithm(Algorithm):
  
  def __init__(self,
               workspaceName: str,
               referencePdbFilename: str,
               predictor: Predictor,
               fitnessFunction: FitnessFunction,
               populationSize: int,
               numIterations: int,
               selection: Selection,
               recombination: Recombination,
               mutation: Mutation,
               replacement: Replacement) -> None:
    super().__init__(workspaceName, referencePdbFilename, predictor)
    self._population_size = populationSize
    self._num_iterations = numIterations
    self._fitness_fn = fitnessFunction
    self._selection = selection
    self._recombination = recombination
    self._mutation = mutation
    self._replacement = replacement
    self._avg_fitnesses = []
  


  @classmethod
  def get_type(cls) -> str:
    return 'GA'



  def get_params_memento(self) -> dict:
    params = super().get_params_memento()
    params['fitnessFunction'] = self._fitness_fn.get_name()
    params['populationSize'] = self._population_size
    params['numIterations'] = self._num_iterations
    params['selection'] = self._selection.get_name()
    params['selectionParams'] = self._selection.get_params_memento()
    params['recombination'] = self._recombination.get_name()
    params['recombinationParams'] = self._recombination.get_params_memento()
    params['mutation'] = self._mutation.get_name()
    params['mutationParams'] = self._mutation.get_params_memento()
    params['replacement'] = self._replacement.get_name()
    params['replacementParams'] = self._replacement.get_params_memento()
    return params

  

  def initialize(self) -> List[Individual]:
    m, n = self._sequence_length, self._population_size
    population = [ Individual.random(m) for _ in range(n) ]
    self.workspace.save_population(0, population)
    self._compute_population_fitness(population)
    return sorted(population)



  def step(self, population: List[Individual]) -> List[Individual]:
    if self.workspace.has_backed_up_children():
      children = self.workspace.restore_children_from_backup()
    else:
      parents = self._selection.apply(population)
      children = self._recombination.apply(parents)
      self._mutation.apply(children)
      self.workspace.backup_children(children)
    for i, child in enumerate(children):
      if child.fitness == None:
        filename = self._get_pdb_filename(child)
        self._fitness_fn.apply(child, 
                              self._predictor, 
                              self._reference_backbone, 
                              filename)
        self.workspace.update_children_backup(i, child)
    self.workspace.delete_children_backup()
    return self._replacement.apply(population, children)



  # TODO: esta funcion, o parte de ella, probablemente se debe mover a la clase 
  # padre
  def run(self, 
          iterationId: int = 0, 
          population: Optional[List[Individual]] = None) -> None:
    if population == None: population = []
    recovering = False
    if len(population) == 0:
      iterationId = 0
      population = self.initialize()
      self.workspace.save_population(iterationId, population)
      stats = Statistics.compute_statistics_and_save(iterationId, 
                                                     population, 
                                                     self.workspace.stats_filename)
      self._avg_fitnesses.append(stats['fitness_mean'])
    else:
      self._load_avg_fitnesses()
      recovering = self._compute_population_fitness(population)
    if recovering:
      population = sorted(population)
      self.workspace.save_population(iterationId, population)
      stats = Statistics.compute_statistics_and_save(iterationId, 
                                                     population, 
                                                     self.workspace.stats_filename)
      self._avg_fitnesses.append(stats['fitness_mean'])
    while True:
      if iterationId == self._num_iterations - 1:
        break
      if population[-1].fitness == self._fitness_fn.upper_bound():
        break
      if iterationId >= 10:
        older_avg = self._avg_fitnesses[iterationId - 10]
        newer_avg = self._avg_fitnesses[iterationId - 1]
        if abs(older_avg - newer_avg) < 0.001:
          break
      population = self.step(population)
      iterationId += 1
      self.workspace.save_population(iterationId, population)
      stats = Statistics.compute_statistics_and_save(iterationId, 
                                                     population, 
                                                     self.workspace.stats_filename)
      self._avg_fitnesses.append(stats['fitness_mean'])
    Statistics.plot_fitness_over_iterations(self.workspace.stats_filename, 
                                            self.workspace.graph_filename,
                                            self.workspace.name)



  def _compute_population_fitness(self, population: List[Individual]) -> bool:
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



  def _load_avg_fitnesses(self) -> None:
    temp = Statistics.load_statistics_from_csv_file(
      self.workspace.stats_filename)
    self._avg_fitnesses = [ stats['fitness_mean'] for stats in temp ]
