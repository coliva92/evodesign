from ..Algorithm import Algorithm
from typing import List, Optional, Tuple
from ..Fitness import FitnessFunction
from ..Prediction import Predictor
from .Selection import Selection
from .Recombination import Recombination
from .Mutation import Mutation
from .Replacement import GA_Replacement_Generational
from ..Individual import Individual
from ..Population import Population
from ..Statistics import Statistics





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
               populationFilenames: Optional[List[str]] = None,
               seed: Optional[int] = None
               ) -> None:
    super().__init__(workspaceName, 
                     targetPdbFilename, 
                     predictor, 
                     fitnessFunction,
                     populationFilenames,
                     seed)
    self._population_size = populationSize
    self._num_iterations = numIterations
    self._selection = selection
    self._recombination = recombination
    self._mutation = mutation
    self._replacement = GA_Replacement_Generational()
    self._terminators = []
  


  def _get_params_json(self) -> dict:
    params = super()._get_params_json()
    params['populationSize'] = self._population_size
    params['numIterations'] = self._num_iterations
    params['selection'] = self._selection.get_name()
    params['selectionParams'] = self._selection.as_json()
    params['recombination'] = self._recombination.get_name()
    params['recombinationParams'] = self._recombination.as_json()
    params['mutation'] = self._mutation.get_name()
    params['mutationParams'] = self._mutation.as_json()
    return params
  


  def _evolutionary_step(self, population: Population) -> List[Individual]:
    new_individuals = []
    while len(new_individuals) < len(population):
      parents = self._selection(population)
      children = self._recombination(parents)
      self._mutation(children)
      new_individuals += children
    return new_individuals



  def next_population(self, population: Population) -> Population:
    children = self.workspace.restore_children_from_backup()
    if not children:
      children = Population(self._evolutionary_step(population))
      self.workspace.backup_children(children)
    try:
      children.update_fitness(self._fitness_fn, 
                              self._predictor, 
                              self._reference_backbone, 
                              self.workspace.pdbs_folder)
    except RuntimeError as e:
      self.workspace.backup_children(children)
      raise e
    self.workspace.delete_children_backup()
    return self._replacement(population, children.individuals)



  def __call__(self, population: Optional[Population] = None) -> None:
    if population is None: population = Population()
    if not population:
      population = Population.new_random(1, 
                                         self._population_size, 
                                         self._sequence_length)
      self.workspace.save_population(population)
    try:
      is_recovering = population.update_fitness(self._fitness_fn,
                                                self._predictor,
                                                self._reference_backbone,
                                                self.workspace.pdbs_folder)
    except RuntimeError as e:
      self.workspace.save_population(population)
      raise e
    self.best_solution = population[-1]
    if is_recovering:
      stats = Statistics.new_from_population(population)
      self.workspace.save_population(population)
      self.workspace.save_statistics(stats, self.best_solution)
    while True:
      if population.iterationId == self._num_iterations:
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
    else:
      population.individuals = population[1:] + [ self.best_solution ]
