from ..FitnessFunction import FitnessFunction
from ...Metrics.Experimental.Cyclization import Cyclization as CycMetric
from ...Metrics.Rmsd import Rmsd
import numpy as np
from scipy.stats import norm





class RmsdCyclization(FitnessFunction):

  _MEAN = 1.3248119
  _STANDARD_DEV = 0.10498072



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
    super().__init__(upperBound, [ Rmsd(), CycMetric() ])
    self._rmsd_weight = rmsdWeight
    self._cyc_weight = cycWeight
    self._weights = np.array([ rmsdWeight, cycWeight ])



  def compute_fitness(self, **kwargs) -> float:
    c = self._compute_probability(kwargs['cyclization'])
    r = self._normalize_rmsd(kwargs['rmsd'])
    return np.average(np.array([ r, c ]), weights=self._weights)
  


  def _compute_probability(self, cyclization: float) -> float:
    z_score = (cyclization - self._MEAN) / self._STANDARD_DEV
    return norm.cdf(z_score)
  


  def _normalize_rmsd(self, rmsd: float) -> float:
    return 1.0 / (1.0 + rmsd)
  