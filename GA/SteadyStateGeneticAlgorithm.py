from .SimpleGeneticAlgorithm import SimpleGeneticAlgorithm
from typing import List, Optional
from ..Fitness import FitnessFunction
from ..Prediction import Predictor
from .Selection import Selection
from .Recombination import Recombination
from .Mutation import Mutation
from .Replacement import Replacement
from .Termination import MeanFitnessStagnation
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
    self._terminators = [ MeanFitnessStagnation(self.workspace.stats_filename) ]



  def _get_params_dict(self) -> dict:
    params = super()._get_params_dict()
    params['replacement'] = self._replacement.get_name()
    return params



  def _evolutionary_step(self, 
                         population: List[Individual]
                         ) -> List[Individual]:
    parents = self._selection.apply(population)
    children = self._recombination(parents)
    self._mutation(children)
    return children
