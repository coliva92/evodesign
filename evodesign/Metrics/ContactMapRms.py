from .Metric import Metric
from typing import List, Optional, Dict
from Bio.PDB.Atom import Atom
import itertools
import numpy as np
import pandas as pd
from ..Context import Context





class ContactMapRms(Metric):
  
  def __init__(self, column: Optional[str] = None) -> None:
    super().__init__(column)
    self._ref_map = None
  


  def _compute_values(self, 
                      backbone: List[Atom],
                      data: pd.Series,
                      context: Context
                      ) -> pd.Series:
    # TODO: guardar los mapas de contacto en un archivo
    if self._ref_map is None:
      self._ref_map = self._contact_map(context.ref_backbone)
    model_map = self._contact_map(backbone)
    difference = self._ref_map - model_map
    rms = np.sqrt(np.mean(difference**2))
    data[self.column_name()] = rms
    return data
  


  def _contact_map(self, backbone: List[Atom]):
    map_values = np.array([ 
      a - b 
      for (a, b) in itertools.combinations(backbone, 2) 
    ])
    return map_values
