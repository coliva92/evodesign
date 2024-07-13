from ..Metric import Metric
from ..Cyclization import Cyclization as CyclizationMetric
from typing import Optional, List, Dict
from ...Context import Context
from Bio.PDB.Atom import Atom
import pandas as pd





class Cyclization(Metric):

  _MEAN = 1.3248119
  _STANDARD_DEVIATION = 0.10498072
  


  def _params(self) -> dict:
    params = super()._params()
    params['cyclization_metric'] = self._cyclization_metric.settings()
    return params
  


  def __init__(self, 
               cyclization_metric: Optional[CyclizationMetric] = None,
               column: Optional[str] = None
               ) -> None:
    super().__init__(column)
    if cyclization_metric is None:
      cyclization_metric = CyclizationMetric()
    self._cyclization_metric = cyclization_metric
  


  def _compute_values(self, 
                      backbone: List[Atom],
                      data: pd.Series,
                      context: Context
                      ) -> pd.Series:
    # compute the cyclization if not already computed
    data = self._cyclization_metric(backbone, data, context)
    z_score =  self._z_score(data[self._cyclization_metric.column_name()])
    data[self.column_name()] = z_score
    return data

    
  
  def _z_score(self, x: float) -> float:
    z = abs(x - self._MEAN) / self._STANDARD_DEVIATION
    return 1.0 / (1.0 + z)
