from ..FitnessFunction import FitnessFunction
from ...Metrics.Experimental.Cyclization import Cyclization as CycMetric
from ...Metrics.Rmsd import Rmsd
import numpy as np





class Cyclization(FitnessFunction):

  @classmethod
  def _class_name(cls) -> str:
    return 'Fitness.Experimental.Cyclization'
  


  @classmethod
  def column_name(cls) -> str:
    return 'fitness_cyclization'
  


  def _params(self) -> dict:
    params = super()._params()
    params['rmsdBound'] = self._rmsd_bound
    params['cyclizationBound'] = self._cyc_bound
    params['rmsdWeight'] = self._rmsd_weight
    params['cycWeight'] = self._cyc_weight
    return params
  


  def __init__(self, 
               upperBound: float = 1.0,
               rmsdBound: float = 2.0,
               cyclizationBound: float = 1.32,
               rmsdWeight: float = 1.0,
               cycWeight: float = 1.0
               ) -> None:
    super().__init__(upperBound, [ CycMetric(), Rmsd() ])
    self._rmsd_bound = rmsdBound
    self._cyc_bound = cyclizationBound
    self._rmsd_weight = rmsdWeight
    self._cyc_weight = cycWeight
    self._weights = np.array([ rmsdWeight, cycWeight ])



  def compute_fitness(self, **kwargs) -> float:
    r = kwargs['rmsd']
    c = kwargs['cyclization']
    r_min = self._rmsd_bound
    c_min = self._cyc_bound
    a = 1.0 / (1.0 + (r - r_min))
    b = 1.0 / (1.0 + (c - c_min))
    return np.average(np.array([ a, b ]), weights=self._weights)
  