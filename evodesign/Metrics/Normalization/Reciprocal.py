from ..Metric import Metric





class Reciprocal(Metric):

  @classmethod
  def _class_name(cls) -> str:
    return 'Metrics.Normalization.Reciprocal'
  


  def column_name(self) -> str:
    return f'reciprocal_norm_{self._metric.column_name()}'
  


  def __init__(self, metric: Metric) -> None:
    super().__init__()
    self._metric = metric
  


  def compute_values(self, **kwargs) -> float:
    otherMetrics = kwargs['otherMetrics']
    value = self._metric(**kwargs)
    otherMetrics[self._metric.column_name()] = value
    return self._normalize(value)
    
  

  def _normalize(self, x: float) -> float:
    return 1.0 / (1.0 + x)
  