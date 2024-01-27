from typing import Dict, Optional
from ..FitnessFunction import FitnessFunction
from ...Metrics.Experimental.Cyclization import Cyclization as CycMetric





class Cyclization(FitnessFunction):

  @classmethod
  def _class_name(cls) -> str:
    return 'Fitness.Cyclization'
  


  @classmethod
  def column_name(cls) -> str:
    return 'Fitness_Cyclization'
  


  @classmethod
  def upper_bound(cls) -> float:
    return -1.4
  


  def __init__(self) -> None:
    super().__init__([ CycMetric() ])



  def compute_fitness(self, 
                      metricValues: Dict[str, float],
                      sequence: Optional[str] = None
                      ) -> float:
    return metricValues['Cyclization']
  