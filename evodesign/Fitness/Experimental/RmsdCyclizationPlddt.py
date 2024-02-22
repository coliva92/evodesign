from ..FitnessFunction import FitnessFunction
from ...Metrics.Experimental.Cyclization import Cyclization
from ...Metrics.Rmsd import Rmsd
import Utils
import numpy as np





class RmsdCyclizationPlddt(FitnessFunction):

  @classmethod
  def _class_name(cls) -> str:
    return 'Fitness.Experimental.RmsdCyclizationPlddt'
  


  @classmethod
  def column_name(cls) -> str:
    return 'fitness_rmsd_cyclization_plddt'
  


  def _params(self) -> dict:
    params = super()._params()
    params['rmsdWeight'] = self._rmsd_weight
    params['cycWeight'] = self._cyc_weight
    params['plddtWeight'] = self._plddt_weight
    return params
  


  def __init__(self, 
               upperBound: float = 1.0,
               rmsdWeight: float = 1.0,
               cycWeight: float = 1.0,
               plddtWeight: float = 1.0
               ) -> None:
    super().__init__(upperBound, [ Rmsd(), Cyclization() ])
    self._rmsd_weight = rmsdWeight
    self._cyc_weight = cycWeight
    self._plddt_weight = plddtWeight
    self._weights = np.array([ rmsdWeight, cycWeight, plddtWeight ])



  def compute_fitness(self, **kwargs) -> float:
    c = Utils.cyclization_probability(kwargs['cyclization'])
    r = Utils.normalize_rmsd(kwargs['rmsd'])
    return np.average(np.array([ r, c, kwargs['plddt'] ]), 
                      weights=self._weights)
