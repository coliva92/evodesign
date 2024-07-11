from ..Metric import Metric
from typing import Optional





class Reciprocal(Metric):
  
  def _params(self) -> dict:
    params = super()._params()
    params['metric'] = self._metric.settings()
    return params
  
  

  def __init__(self, 
               metric: Metric,
               column: Optional[str] = None
               ) -> None:
    super().__init__(column)
    self._metric = metric
  


  def compute_value(self, **kwargs) -> float:
    otherMetrics = kwargs['otherMetrics']
    value = self._metric(**kwargs)
    otherMetrics[self._metric.column_name()] = value
    return self._normalize(value)
    
  

  def _normalize(self, x: float) -> float:
    return 1.0 / (1.0 + x)
  