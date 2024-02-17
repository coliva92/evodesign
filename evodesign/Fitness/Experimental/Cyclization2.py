from .Cyclization import Cyclization
import numpy as np





class Cyclization2(Cyclization):

  @classmethod
  def _class_name(cls) -> str:
    return 'Fitness.Experimental.Cyclization2'
  


  @classmethod
  def column_name(cls) -> str:
    return 'fitness_cyclization2'



  def compute_fitness(self, **kwargs) -> float:
    r = kwargs['rmsd']
    c = kwargs['cyclization']
    r_min = self._rmsd_bound
    c_min = self._cyc_bound
    return np.average(np.array([ r_min / (1.0 + r), 1.0 - c + c_min ]), 
                      weights=self._weights)
  