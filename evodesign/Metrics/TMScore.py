from .Metric import Metric
from .Rmsd import Rmsd
import numpy as np
from typing import Optional, Dict, List
from ..Context import Context
from Bio.PDB.Atom import Atom
import pandas as pd





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
    if rmsd_metric is None:
      rmsd_metric = Rmsd()
    self._rmsd_metric = rmsd_metric
  


  def _compute_values(self, 
                     backbone: List[Atom],
                     data: pd.Series,
                     context: Context
                     ) -> pd.Series:
    # superimpose the structures are not already superimposed
    data = self._rmsd_metric(backbone, data, context)
    d0 = self._normalizing_constant(len(context.ref_backbone))
    distances = np.array([ a - b for a, b in zip(backbone, context.ref_backbone) ])
    tmscore = np.mean(1 / (1 + (distances / d0) ** 2))
    data[self.column_name()] = tmscore
    return data
    
  

  def _normalizing_constant(self, n: int) -> float:
    return 1.24 * np.cbrt(n - 15) - 1.8
