from ..FitnessFunction import FitnessFunction
from ...Metrics.Experimental.Cyclization import Cyclization as CycMetric





class Cyclization(FitnessFunction):

  @classmethod
  def _class_name(cls) -> str:
    return 'Fitness.Experimental.Cyclization'
  


  @classmethod
  def column_name(cls) -> str:
    return 'fitness_cyclization'
  


  def __init__(self, upperBound: float = -1.4) -> None:
    super().__init__(upperBound, [ CycMetric() ])



  def compute_fitness(self, **kwargs) -> float:
    return -kwargs['cyclization']
  