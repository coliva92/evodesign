from .SimpleGeneticAlgorithm import SimpleGeneticAlgorithm
from ..Fitness import FitnessFunction
from ..Prediction import Predictor
from .Selection import Selection
from .Recombination import Recombination
from .Mutation import Mutation
from .Replacement import Replacement
from .Termination import StagnatingMeanFitness
from ..Population import Population





class SteadyStateGeneticAlgorithm(SimpleGeneticAlgorithm):

  @classmethod
  def get_name(cls) -> str:
    return 'GA_Steady'
  

  
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
               replacement: Replacement
               ) -> None:
    super().__init__(workspaceRoot, 
                     targetPdbFilename, 
                     predictor, 
                     fitnessFunction, 
                     populationSize, 
                     numIterations, 
                     selection, 
                     recombination, 
                     mutation)
    self._replacement = replacement
    self._terminators = [ StagnatingMeanFitness(self.workspace.stats_filename) ]



  def _params_json(self) -> dict:
    params = super()._params_json()
    params['replacement'] = self._replacement.get_name()
    return params



  def _evolutionary_step(self, population: Population) -> Population:
    parents = self._selection(population)
    children = self._recombination(parents)
    self._mutation(children)
    return children
