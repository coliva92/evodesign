from .SteadyStateGeneticAlgorithm import SteadyStateGeneticAlgorithm
from typing import List
from evodesign.Fitness import FitnessFunction
from evodesign.GA.Mutation import Mutation
from evodesign.GA.Recombination import Recombination
from evodesign.GA.Selection import Selection
from evodesign.Prediction import Predictor
from ..Population import Individual





class SimpleGeneticAlgorithm(SteadyStateGeneticAlgorithm):
  
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
    super().__init__(workspaceName, 
                     targetPdbFilename, 
                     predictor, 
                     fitnessFunction, 
                     populationSize, 
                     numIterations, 
                     selection, 
                     recombination, 
                     mutation, 
                     None)
  


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
