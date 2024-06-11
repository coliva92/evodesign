from ..Metric import Metric
from ..Cyclization import Cyclization





class Cyclization(Metric):

  _MEAN = 1.3248119
  _STANDARD_DEVIATION = 0.10498072



  @classmethod
  def _class_name(cls) -> str:
    return 'Metrics.Normalization.Cyclization'
  


  def column_name(self) -> str:
    return 'norm_cyclization'
  


  def _params(self) -> dict:
    return {
      'cycCalc': self._cyclization.settings()
    }
  


  def __init__(self, cycCalc: Cyclization) -> None:
    super().__init__()
    self._cyclization = cycCalc
  


  def compute_value(self, **kwargs) -> float:
    otherMetrics = kwargs['otherMetrics']
    value = self._cyclization(**kwargs)
    otherMetrics[self._cyclization.column_name()] = value
    return self._z_score(value)

    
  
  def _z_score(self, x: float) -> float:
    z = abs(x - self._MEAN) / self._STANDARD_DEVIATION
    return z
