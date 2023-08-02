from .SimpleGeneticAlgorithm import SimpleGeneticAlgorithm
from typing import List, Optional
from ..Fitness import FitnessFunction
from ..Prediction import Predictor
from .Selection import Selection
from .Recombination import Recombination
from .Mutation import Mutation
from .Replacement import Replacement
from ..Individual import Individual





class SteadyStateGeneticAlgorithm(SimpleGeneticAlgorithm):

  @classmethod
  def get_name(cls) -> str:
    return 'GA_Steady'
  

  
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
               replacement: Replacement,
               populationFilenames: Optional[List[str]] = None
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
                     populationFilenames)
    self._replacement = replacement
    self._avg_fitnesses = []



  def _get_params_memento(self) -> dict:
    params = super()._get_params_memento()
    params['replacement'] = self._replacement.get_name()
    params['replacementParams'] = self._replacement.get_params_memento()
    return params



  def _evolutionary_step(self, 
                         population: List[Individual]
                         ) -> List[Individual]:
    parents = self._selection.apply(population)
    children = self._recombination.apply(parents)
    self._mutation.apply(children)
    return children



  def _additional_termination_conditions(self, iterationId: int) -> bool:
    if iterationId >= 10:
      older_avg = self._statistics[iterationId - 10].fitness_mean
      newer_avg = self._statistics[iterationId].fitness_mean
      return newer_avg - older_avg < 0.0001
