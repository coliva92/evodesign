from .Metric import Metric
from .Rmsd import Rmsd
import numpy as np





class TMScore(Metric):

  @classmethod
  def _class_name(cls) -> str:
    return 'Metrics.TMScore'
  


  def column_name(self) -> str:
    return 'tm_score'
  


  def _params(self) -> dict:
    return {
      'rmsdCalculator': self._rmsd_calc.settings()
    }
  


  def __init__(self,
               rmsdCalculator: Rmsd
               ) -> None:
    super().__init__()
    self._rmsd_calc = rmsdCalculator
  


  def compute_value(self, **kwargs) -> float:
    c = self._rmsd_calc.column_name()
    # check if the structures were already superimposed
    if c not in kwargs['otherMetrics']:
      kwargs['otherMetrics'][c] = self._rmsd_calc.compute_value(**kwargs)
    model, reference = kwargs['model'], kwargs['reference']
    d0 = self._normalizing_constant(len(reference))
    distances = np.array([ a - b for a, b in zip(model, reference) ])
    return np.mean(1 / (1 + (distances / d0)**2))
    
  

  def _normalizing_constant(self, n: int) -> float:
    return 1.24 * np.cbrt(n - 15) - 1.8
