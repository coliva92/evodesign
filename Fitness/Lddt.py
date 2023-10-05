from .FitnessFunction import FitnessFunction
from typing import Dict, List
from ..Metrics import Lddt





class Lddt(FitnessFunction):

  @classmethod
  def name(cls) -> str:
    return 'Fitness_Lddt'
  


  @classmethod
  def upper_bound(cls) -> float:
    return 0.95
  


  def __init__(self,
               cutoffs: List[float] = [ 0.5, 1.0, 2.0, 4.0 ],
               inclusionRadius: float = 1
               ) -> None:
    metrics = {
      'lddt': Lddt(cutoffs, inclusionRadius)
    }
    super().__init__(metrics)
  


  def params_json(self) -> dict:
    return {
      'cutoffs': self._metric_calculators['lddt']._cutoffs,
      'inclusionRadius': self._metric_calculators['lddt']._radius
    }
  


  def compute_fitness(self, metrics: Dict[str, float]) -> float:
    return metrics['lddt']
