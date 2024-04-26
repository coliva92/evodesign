from ..FitnessFunction import FitnessFunction
from typing import List
from ...Metrics.Experimental.Cyclization import Cyclization
from ...Metrics.Rmsd import Rmsd
from ...Metrics.Gdt import Gdt
from ...Metrics.Experimental.ContactMapRMS import ContactMapRMS
import evodesign.Fitness.Experimental.Utils as Utils
import numpy as np





class GdtContactMapCyclizationPlddt(FitnessFunction):

  @classmethod
  def _class_name(cls) -> str:
    return 'Fitness.Experimental.GdtContactMapCyclizationPlddt'
  


  @classmethod
  def column_name(cls) -> str:
    return 'fitness_gdt_cm_cyclization_plddt'
  


  def _params(self) -> dict:
    params = super()._params()
    params['cutoffs'] = self._cutoffs
    params['gdtWeight'] = self._gdt_weight
    params['cmWeight'] = self._contact_map_weight
    params['cycWeight'] = self._cyc_weight
    params['plddtWeight'] = self._plddt_weight
    return params
  


  def __init__(self, 
               upperBound: float = 1.0,
               cutoffs: List[float] = [ 1.0, 2.0, 4.0, 8.0 ],
               gdtWeight: float = 1.0,
               cmWeight: float = 1.0,
               cycWeight: float = 1.0,
               plddtWeight: float = 1.0
               ) -> None:
    metrics = [ Rmsd(), Gdt(cutoffs), ContactMapRMS(), Cyclization() ]
    super().__init__(upperBound, metrics)
    self._cutoffs = cutoffs
    self._gdt_weight = gdtWeight
    self._cyc_weight = cycWeight
    self._contact_map_weight = cmWeight
    self._plddt_weight = plddtWeight
    self._weights = np.array([ gdtWeight, cmWeight, cycWeight, plddtWeight ])

  


  def compute_fitness(self, **kwargs) -> float:
    cyc = Utils.normalize(abs(Utils.cyclization_z_score(kwargs['cyclization'])))
    plddt = Utils.normalize_plddt(kwargs['plddt'])
    cm_rms = Utils.normalize(kwargs['contact_maps_rms'])
    terms = np.array([ kwargs['gdt'], cm_rms, cyc, plddt ])
    return np.average(terms, weights=self._weights)
