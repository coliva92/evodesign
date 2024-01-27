from ..FitnessFunction import FitnessFunction
from typing import Dict, Optional
from ...Metrics.Experimental.SideChainPacking import SideChainPacking \
  as PackingMetric
import math





class SideChainPacking(FitnessFunction):

  @classmethod
  def _class_name(cls) -> str:
    return 'Fitness.SideChainPacking'
  


  @classmethod
  def column_name(cls) -> str:
    return 'Side_Chain_Packing'
  


  def _params(self) -> dict:
    return self._metrics[0].params()
  


  @classmethod
  def upper_bound(cls) -> float:
    return math.inf
  


  def __init__(self,
               scwrlExecutablePath: str = './scwrl4/Scwrl4') -> None:
    super().__init__([ PackingMetric(scwrlExecutablePath) ])
  


  def compute_fitness(self, 
                      metricValues: Dict[str, float], 
                      sequence: Optional[str] = None
                      ) -> float:
    return -metricValues['Side_Chain_Packing']
  