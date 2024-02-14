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
  


  def __init__(self, 
               upperBound: float = 1.0,
               rmsdBound: float = 2.0,
               cyclizationBound: float = 1.5) -> None:
    super().__init__(upperBound, [ CycMetric(), Rmsd() ])
    self._rmsd_bound = rmsdBound
    self._cyc_bound = cyclizationBound



  def compute_fitness(self, **kwargs) -> float:
    r = kwargs['rmsd']
    c = kwargs['cyclization']
    r_min = self._rmsd_bound
    c_min = self._cyc_bound
    return np.array([ r_min / r, c_min / c ])
  