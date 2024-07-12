from ..Metric import Metric
from typing import Optional, List, Dict
from ...Context import Context
from Bio.PDB.Atom import Atom
import pandas as pd





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
  


  def _compute_values(self, 
                      backbone: List[Atom],
                      data: pd.Series,
                      context: Context
                      ) -> pd.Series:
    # compute the metric if not already computed
    data = self._metric(backbone, data, context)
    value = data[self._metric.column_name()]
    data[self.column_name()] = self._normalize(value)
    return data
    
  

  def _normalize(self, x: float) -> float:
    return 1.0 / (1.0 + x)
  