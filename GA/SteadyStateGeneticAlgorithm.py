from .SimpleGeneticAlgorithm import SimpleGeneticAlgorithm
from typing import List, Optional
from ..Fitness import FitnessFunction
from ..Prediction import Predictor
from .Selection import Selection
from .Recombination import Recombination
from .Mutation import Mutation
from .Replacement import Replacement
from ..Population import Individual
import evodesign.Statistics as Statistics





class SteadyStateGeneticAlgorithm(SimpleGeneticAlgorithm):

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
               replacement: Replacement
               ) -> None:
    super().__init__(workspaceName, 
                     targetPdbFilename, 
                     predictor, 
                     fitnessFunction, 
                     populationSize, 
                     numIterations, 
                     selection, 
                     recombination, 
                     mutation)
    self._replacement = replacement
    self._avg_fitnesses = []



  @classmethod
  def get_type(cls) -> str:
    return 'GA_Steady'



  def _get_params_memento(self) -> dict:
    params = super()._get_params_memento()
    params['replacement'] = self._replacement.get_name()
    params['replacementParams'] = self._replacement.get_params_memento()
    return params



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
    is_recovering = False
    if len(population) == 0:
      iterationId = 0
      population = self.initialize()
      stats = Statistics.compute_statistics(population)
      self._avg_fitnesses.append(stats['fitness_mean'])
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
      self._avg_fitnesses.append(stats['fitness_mean'])
      self.workspace.save_population(iterationId, population)
      Statistics.save_statistics_to_csv_file(stats, 
                                             iterationId, 
                                             self.workspace.stats_filename)
    while True:
      if iterationId == self._num_iterations - 1:
        break
      if population[-1].fitness == self._fitness_fn.upper_bound():
        break
      if iterationId >= 10 and self.get_type() == 'GA_Steady':
        older_avg = self._avg_fitnesses[iterationId - 10]
        newer_avg = self._avg_fitnesses[iterationId]
        if newer_avg - older_avg < 0.0001:
          break
      population = self.step(population)
      iterationId += 1
      stats = Statistics.compute_statistics(population)
      self._avg_fitnesses.append(stats['fitness_mean'])
      self.workspace.save_population(iterationId, population)
      Statistics.save_statistics_to_csv_file(stats, 
                                             iterationId, 
                                             self.workspace.stats_filename)
    Statistics.plot_fitness_over_iterations(self.workspace.stats_filename, 
                                            self.workspace.graph_filename,
                                            self.workspace.name)



  def _recovering_statistics(self) -> None:
    temp = Statistics.load_statistics_from_csv_file(
      self.workspace.stats_filename)
    self._avg_fitnesses = [ stats['fitness_mean'] for stats in temp ]
  