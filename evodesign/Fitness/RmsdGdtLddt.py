from . import FitnessFunction
from typing import Dict, List, Optional
from ..Metrics import Rmsd, Gdt, Lddt





class RmsdGdtLddt(FitnessFunction):

  @classmethod
  def name(cls) -> str:
    return 'Fitness_RmsdGdtLddt'
  


  @classmethod
  def upper_bound(cls) -> float:
    return 0.95
  


  def __init__(self,
               gdtCutoffs: List[float] = [ 1.0, 2.0, 4.0, 8.0 ],
               lddtCutoffs: List[float] = [ 0.5, 1.0, 2.0, 4.0 ],
               inclusionRadius: float = 15
               ) -> None:
    metrics = {
      'rmsd': Rmsd(),
      'gdt': Gdt(gdtCutoffs),
      'lddt': Lddt(lddtCutoffs, inclusionRadius)
    }
    super().__init__(metrics)
  


  def params_as_dict(self) -> dict:
    return {
      'gdtCutoffs': self._metric_calculators['gdt']._cutoffs,
      'lddtCutoffs': self._metric_calculators['lddt']._cutoffs
    }
  


  def compute_fitness(self, 
                      metrics: Dict[str, float],
                      _: Optional[str]) -> float:
    return ((1.95 / metrics['rmsd']) + metrics['gdt'] + metrics['lddt']) / 3
