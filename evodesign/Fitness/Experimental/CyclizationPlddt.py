from ..FitnessFunction import FitnessFunction
from ...Metrics.Experimental.Cyclization import Cyclization
from ...Metrics.Rmsd import Rmsd
import evodesign.Fitness.Experimental.Utils as Utils
import numpy as np





class CyclizationPlddt(FitnessFunction):

  @classmethod
  def _class_name(cls) -> str:
    return 'Fitness.Experimental.CyclizationPlddt'
  


  @classmethod
  def column_name(cls) -> str:
    return 'fitness_cyclization_plddt'
  


  def _params(self) -> dict:
    params = super()._params()
    params['cycWeight'] = self._cyc_weight
    params['plddtWeight'] = self._plddt_weight
    return params
  


  def __init__(self, 
               upperBound: float = 1.0,
               cycWeight: float = 1.0,
               plddtWeight: float = 1.0
               ) -> None:
    super().__init__(upperBound, [ Rmsd(), Cyclization() ])
    self._cyc_weight = cycWeight
    self._plddt_weight = plddtWeight
    self._weights = np.array([ cycWeight, plddtWeight ])



  def compute_fitness(self, **kwargs) -> float:
    c = Utils.normalize(abs(Utils.cyclization_z_score(kwargs['cyclization'])))
    p = Utils.normalize_plddt(kwargs['plddt'])
    return np.average(np.array([ c, p ]), weights=self._weights)
