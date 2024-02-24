from ..FitnessFunction import FitnessFunction
from ...Metrics.Experimental.Cyclization import Cyclization
from ...Metrics.Rmsd import Rmsd
import evodesign.Fitness.Experimental.Utils as Utils
import numpy as np





class RmsdCyclization(FitnessFunction):

  @classmethod
  def _class_name(cls) -> str:
    return 'Fitness.Experimental.RmsdCyclization'
  


  @classmethod
  def column_name(cls) -> str:
    return 'fitness_rmsd_cyclization'
  


  def _params(self) -> dict:
    params = super()._params()
    params['rmsdWeight'] = self._rmsd_weight
    params['cycWeight'] = self._cyc_weight
    return params
  


  def __init__(self, 
               upperBound: float = 1.0,
               rmsdWeight: float = 1.0,
               cycWeight: float = 1.0
               ) -> None:
    super().__init__(upperBound, [ Rmsd(), Cyclization() ])
    self._rmsd_weight = rmsdWeight
    self._cyc_weight = cycWeight
    self._weights = np.array([ rmsdWeight, cycWeight ])



  def compute_fitness(self, **kwargs) -> float:
    c = Utils.normalize_cyclization(kwargs['cyclization'])
    r = Utils.normalize_rmsd(kwargs['rmsd'])
    return np.average(np.array([ r, c ]), weights=self._weights)
