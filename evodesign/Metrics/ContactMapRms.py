from .Metric import Metric
from typing import List
from Bio.PDB.Atom import Atom
import itertools
import numpy as np





class ContactMapRms(Metric):

  _reference_map = None
  


  def column_name(self) -> str:
    return 'contact_map_rms'
  


  def compute_value(self, **kwargs) -> float:
    # TODO: guardar los mapas de contacto en un archivo
    model, reference = kwargs['model'], kwargs['reference']
    if self._reference_map is None:
      self._reference_map = self._contact_map(reference)
    model_map = self._contact_map(model)
    difference = self._reference_map - model_map
    rms = np.sqrt(np.mean(difference**2))
    return rms
  


  def _contact_map(self, backbone: List[Atom]):
    map_values = np.array([ 
      a - b 
      for (a, b) in itertools.combinations(backbone, 2) 
    ])
    return map_values
