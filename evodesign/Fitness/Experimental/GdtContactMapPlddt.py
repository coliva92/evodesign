from typing import List
from ..FitnessFunction import FitnessFunction
from Metrics.Rmsd import Rmsd
from Metrics.Gdt import Gdt
from Metrics.Experimental.ContactMapRMS import ContactMapRMS
import evodesign.Fitness.Experimental.Utils as Utils
import numpy as np





class GdtContactMapPlddt(FitnessFunction):

  @classmethod
  def _class_name(cls) -> str:
    return 'Fitness.Experimental.GdtContactMap'
  


  @classmethod
  def column_name(cls) -> str:
    return 'fitness_gdt_cmap'
  


  def _params(self) -> dict:
    params = super()._params()
    params['cutoffs'] = self._cutoffs
    params['gdtWeight'] = self._gdt_weight
    params['contactMapWeight'] = self._contact_map_weight
    params['plddtWeight'] = self._plddt_weight
    return params
  


  def __init__(self, 
               upperBound: float = 1.0,
               cutoffs: List[float] = [ 1.0, 2.0, 4.0, 8.0 ],
               gdtWeight: float = 1.0,
               contactMapWeight: float = 1.0,
               plddtWeight: float = 1.0
               ) -> None:
    super().__init__(upperBound, [ Rmsd(), Gdt(cutoffs), ContactMapRMS() ])
    self._cutoffs = cutoffs
    self._gdt_weight = gdtWeight
    self._contact_map_weight = contactMapWeight
    self._plddt_weight = plddtWeight
    self._weights = np.array([ gdtWeight, contactMapWeight, plddtWeight ])

  


  def compute_fitness(self, **kwargs) -> float:
    plddt = Utils.normalize_plddt(kwargs['plddt'])
    terms = np.array([ kwargs['gdt'], kwargs['contact_maps_rms'], plddt ])
    return np.average(terms, weights=self._weights)
