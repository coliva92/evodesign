from .Metric import Metric
from .Rmsd import Rmsd
import numpy as np
from typing import Optional





class TMScore(Metric):

  def _params(self) -> dict:
    params = super()._params()
    params['rmsd_metric'] = self._rmsd_metric.settings()
    return params



  def __init__(self, 
               rmsd_metric: Optional[Rmsd] = None,
               column: Optional[str] = None
               ) -> None:
    super().__init__(column)
    self._rmsd_metric = rmsd_metric
    if self._rmsd_metric is None:
      self._rmsd_metric = Rmsd()
  


  def compute_value(self, **kwargs) -> float:
    c = self._rmsd_metric.column_name()
    # check if the structures were already superimposed
    if c not in kwargs['otherMetrics']:
      kwargs['otherMetrics'][c] = self._rmsd_metric.compute_value(**kwargs)
    model, reference = kwargs['model'], kwargs['reference']
    d0 = self._normalizing_constant(len(reference))
    distances = np.array([ a - b for a, b in zip(model, reference) ])
    return np.mean(1 / (1 + (distances / d0)**2))
    
  

  def _normalizing_constant(self, n: int) -> float:
    return 1.24 * np.cbrt(n - 15) - 1.8
