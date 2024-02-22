from ..FitnessFunction import FitnessFunction
from ...Metrics.Rmsd import Rmsd
import Utils
import numpy as np





class RmsdPlddt(FitnessFunction):

  @classmethod
  def _class_name(cls) -> str:
    return 'Fitness.Experimental.RmsdPlddt'
  


  @classmethod
  def column_name(cls) -> str:
    return 'fitness_rmsd_plddt'
  


  def _params(self) -> dict:
    params = super()._params()
    params['rmsdWeight'] = self._rmsd_weight
    params['plddtWeight'] = self._plddt_weight
    return params
  


  def __init__(self, 
               upperBound: float = 1.0,
               rmsdWeight: float = 1.0,
               plddtWeight: float = 1.0
               ) -> None:
    super().__init__(upperBound, [ Rmsd() ])
    self._rmsd_weight = rmsdWeight
    self._plddt_weight = plddtWeight
    self._weights = np.array([ rmsdWeight, plddtWeight ])
  


  def compute_fitness(self, **kwargs) -> float:
    r = Utils.normalize_rmsd(kwargs['rmsd'])
    p = Utils.normalize_plddt(kwargs['plddt'])
    return np.average(np.array([ r, p ]), weights=self._weights)
