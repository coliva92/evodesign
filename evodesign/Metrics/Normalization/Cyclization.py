from ..Metric import Metric
from ..Cyclization import Cyclization as CyclizationMetric
from typing import Optional





class Cyclization(Metric):

  _MEAN = 1.3248119
  _STANDARD_DEVIATION = 0.10498072
  


  def _params(self) -> dict:
    params = super()._params()
    params['cyclization_calc'] = self._cyclization_metric.settings()
    return params
  


  def __init__(self, 
               cyclization_metric: Optional[CyclizationMetric] = None,
               column: Optional[str] = None
               ) -> None:
    super().__init__(column)
    if cyclization_metric is None:
      cyclization_metric = CyclizationMetric()
    self._cyclization_metric = cyclization_metric
  


  def compute_value(self, **kwargs) -> float:
    otherMetrics = kwargs['otherMetrics']
    value = self._cyclization_metric(**kwargs)
    otherMetrics[self._cyclization_metric.column_name()] = value
    return self._z_score(value)

    
  
  def _z_score(self, x: float) -> float:
    z = abs(x - self._MEAN) / self._STANDARD_DEVIATION
    return z
