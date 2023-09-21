from .FitnessFunction import FitnessFunction
from typing import Dict, List
from ..Metrics import Rmsd, Gdt, Lddt





class GdtLddt(FitnessFunction):

  @classmethod
  def name(cls) -> str:
    return 'Fitness_GdtLddt'
  


  @classmethod
  def upper_bound(cls) -> float:
    return 1.9
  


  def __init__(self,
               gdtCutoffs: List[float] = [ 1.0, 2.0, 4.0, 8.0 ],
               lddtCutoffs: List[float] = [ 0.5, 1.0, 2.0, 4.0 ]
               ) -> None:
    metrics = {
      'rmsd': Rmsd(),
      'gdt': Gdt(gdtCutoffs),
      'lddt': Lddt(lddtCutoffs)
    }
    super().__init__(metrics)
  


  def params_json(self) -> dict:
    return {
      'gdtCutoffs': self._metric_calculators['gdt']._cutoffs,
      'lddtCutoffs': self._metric_calculators['lddt']._cutoffs
    }
  


  def compute_fitness(self, metrics: Dict[str, float]) -> float:
    return metrics['gdt'] + metrics['lddt']
