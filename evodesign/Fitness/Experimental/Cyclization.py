from ..FitnessFunction import FitnessFunction
from typing import List
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
    params['weights'] = self._weights
    return params
  


  def __init__(self, 
               upperBound: float = 1.0,
               rmsdBound: float = 2.0,
               cyclizationBound: float = 1.32,
               weights: List[float] = [ 1.0, 0.5 ]) -> None:
    super().__init__(upperBound, [ CycMetric(), Rmsd() ])
    self._rmsd_bound = rmsdBound
    self._cyc_bound = cyclizationBound
    self._weights = weights



  def compute_fitness(self, **kwargs) -> float:
    r = kwargs['rmsd']
    c = kwargs['cyclization']
    r_min = self._rmsd_bound
    c_min = self._cyc_bound
    return np.average(np.array([ r_min / (1.0 + r), c_min / (1.0 + c) ]), 
                      weights=self._weights)
  