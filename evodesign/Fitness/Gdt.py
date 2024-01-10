from .FitnessFunction import FitnessFunction
from typing import Dict, List, Optional
from ..Metrics import Rmsd, Gdt as GdtMetric





class Gdt(FitnessFunction):

  @classmethod
  def name(cls) -> str:
    return 'Fitness_Gdt'
  


  @classmethod
  def upper_bound(cls) -> float:
    return 0.95
  


  def __init__(self,
               cutoffs: List[float] = [ 1.0, 2.0, 4.0, 8.0 ],
               ) -> None:
    metrics = {
      'rmsd': Rmsd(),
      'gdt': GdtMetric(cutoffs)
    }
    super().__init__(metrics)
  


  def params_json(self) -> dict:
    return {
      'cutoffs': self._metric_calculators['gdt']._cutoffs
    }
  


  def compute_fitness(self, 
                      metrics: Dict[str, float],
                      _: Optional[str] = None) -> float:
    return metrics['gdt']
