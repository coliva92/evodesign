from .FitnessFunction import FitnessFunction
from typing import Dict, List
from ..Metrics import Rmsd, Gdt as GdtMetric





class GdtRmsdRatio(FitnessFunction):
  
  @classmethod
  def name(cls) -> str:
    return 'Fitness_GdtRmsdRatio'
  


  @classmethod
  def upper_bound(cls) -> float:
    return 95
  


  def __init__(self,
               cutoffs: List[float] = [ 1.0, 2.0, 4.0, 8.0 ],
               carbonAlphaOnly: bool = False
               ) -> None:
    metrics = {
      'rmsd': Rmsd(),
      'gdt': GdtMetric(cutoffs, carbonAlphaOnly)
    }
    super().__init__(metrics)
  


  def params_json(self) -> dict:
    return {
      'cutoffs': self._metric_calculators['gdt']._cutoffs,
      'alphaCarbonOnly': self._metric_calculators['gdt']._carbon_alpha_only
    }
  


  def compute_fitness(self, metrics: Dict[str, float]) -> float:
    return 2 * metrics['gdt'] / metrics['rmsd']
  