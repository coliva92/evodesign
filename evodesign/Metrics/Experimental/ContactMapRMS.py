from ..Metric import Metric
from typing import List
from Bio.PDB.Atom import Atom
import itertools
import numpy as np





class ContactMapRMS(Metric):

  @classmethod
  def column_name(cls) -> str:
    return 'contact_maps_rms'
  


  def __call__(self, **kwargs) -> float:
    model, reference = kwargs['model'], kwargs['reference']
    model_map = self._contact_map(model)
    reference_map = self._contact_map(reference)
    difference = reference_map - model_map
    rms = np.sqrt(np.mean(difference**2))
    return rms
  


  def _contact_map(self, backbone: List[Atom]):
    map_values = np.array([ 
      a - b 
      for (a, b) in itertools.combinations(backbone, 2) 
    ])
    return map_values
