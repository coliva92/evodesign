from ..FitnessFunction import FitnessFunction
from ...Metrics.Experimental.SideChainPacking import SideChainPacking \
  as PackingMetric
import math





class SideChainPacking(FitnessFunction):

  @classmethod
  def _class_name(cls) -> str:
    return 'Fitness.Experimental.SideChainPacking'
  


  @classmethod
  def column_name(cls) -> str:
    return 'fitness_sidechain_packing'
  


  def _params(self) -> dict:
    params = super()._params()
    params['scwrlExecutablePath'] = self._scwrl_executable_path
    return params
  


  def __init__(self,
               upperBound: float = -math.inf,
               scwrlExecutablePath: str = './scwrl4/Scwrl4'
               ) -> None:
    super().__init__(upperBound, [ PackingMetric(scwrlExecutablePath) ])
    self._scwrl_executable_path = scwrlExecutablePath
  


  def compute_fitness(self, **kwargs) -> float:
    return -kwargs['sidechain_packing']
  