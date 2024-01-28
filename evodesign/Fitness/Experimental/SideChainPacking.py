from ..FitnessFunction import FitnessFunction
from ...Metrics.Experimental.SideChainPacking import SideChainPacking \
  as PackingMetric
import math





class SideChainPacking(FitnessFunction):

  @classmethod
  def _class_name(cls) -> str:
    return 'Fitness.SideChainPacking'
  


  @classmethod
  def column_name(cls) -> str:
    return 'fitness_sidechain_packing'
  


  def _params(self) -> dict:
    return self._metrics[0].params()
  


  @classmethod
  def upper_bound(cls) -> float:
    return math.inf
  


  def __init__(self,
               scwrlExecutablePath: str = './scwrl4/Scwrl4') -> None:
    super().__init__([ PackingMetric(scwrlExecutablePath) ])
  


  def compute_fitness(self, **kwargs) -> float:
    return -kwargs['sidechain_packing']
  