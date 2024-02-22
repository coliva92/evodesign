from ..FitnessFunction import FitnessFunction
from typing import List
from ...Metrics.Rmsd import Rmsd
from ...Metrics.Gdt import Gdt
import evodesign.Fitness.Experimental.Utils as Utils
import numpy as np





class GdtPlddt(FitnessFunction):

  @classmethod
  def _class_name(cls) -> str:
    return 'Fitness.Experimental.GdtPlddt'
  


  @classmethod
  def column_name(cls) -> str:
    return 'fitness_gdt_plddt'
  


  def _params(self) -> dict:
    params = super()._params()
    params['cutoffs'] = self._cutoffs
    params['gdtWeight'] = self._gdt_weight
    params['plddtWeight'] = self._plddt_weight
    return params
  


  def __init__(self, 
               upperBound: float = 1.0,
               cutoffs: List[float] = [ 1.0, 2.0, 4.0, 8.0 ],
               gdtWeight: float = 1.0,
               plddtWeight: float = 1.0
               ) -> None:
    super().__init__(upperBound, [ Rmsd(), Gdt(cutoffs) ])
    self._cutoffs = cutoffs
    self._gdt_weight = gdtWeight
    self._plddt_weight = plddtWeight
    self._weights = np.array([ gdtWeight, plddtWeight ])
  


  def compute_fitness(self, **kwargs) -> float:
    p = Utils.normalize_plddt(kwargs['plddt'])
    return np.average(np.array([ kwargs['gdt'], p ]), weights=self._weights)
  