from ..FitnessFunction import FitnessFunction
from ...Metrics.Experimental.Cyclization import Cyclization as CycMetric





class Cyclization(FitnessFunction):

  @classmethod
  def _class_name(cls) -> str:
    return 'Fitness.Cyclization'
  


  @classmethod
  def column_name(cls) -> str:
    return 'fitness_cyclization'
  


  @classmethod
  def upper_bound(cls) -> float:
    return -1.4
  


  def __init__(self) -> None:
    super().__init__([ CycMetric() ])



  def compute_fitness(self, **kwargs) -> float:
    return kwargs['cyclization']
  