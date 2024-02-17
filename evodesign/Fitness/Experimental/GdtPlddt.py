from ..FitnessFunction import FitnessFunction
from typing import List
from ...Metrics.Rmsd import Rmsd
from ...Metrics.Gdt import Gdt





class GdtPlddt(FitnessFunction):

  @classmethod
  def _class_name(cls) -> str:
    return 'Fitness.GdtPlddt'
  


  @classmethod
  def column_name(cls) -> str:
    return 'fitness_gdt_plddt'
  


  def _params(self) -> dict:
    params = super()._params()
    params['cutoffs'] = self._cutoffs
    return params
  


  def __init__(self, 
               upperBound: float = 0.855,
               cutoffs: List[float] = [ 1.0, 2.0, 4.0, 8.0 ]
               ) -> None:
    super().__init__(upperBound, [ Rmsd(), Gdt(cutoffs) ])
    self._cutoffs = cutoffs
  


  def compute_fitness(self, **kwargs) -> float:
    return kwargs['gdt'] * kwargs['plddt']
  