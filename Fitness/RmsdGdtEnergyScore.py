from .FitnessFunction import FitnessFunction
from typing import Dict, List
from ..Metrics import Rmsd, Gdt, EnergyScore
import math





class RmsdGdtEnergyScore(FitnessFunction):

  @classmethod
  def name(cls) -> str:
    return 'Fitness_RmsdGdtEnergyScore'
  


  @classmethod
  def upper_bound(cls) -> float:
    return 3.5
  


  def __init__(self,
               gdtCutoffs: List[float] = [ 1.0, 2.0, 4.0, 8.0 ]
               ) -> None:
    metrics = {
      'rmsd': Rmsd(),
      'gdt': Gdt(gdtCutoffs),
      'energyScore': EnergyScore()
    }
    super().__init__(metrics)
  


  def params_json(self) -> dict:
    return {
      'gdtCutoffs': self._metric_calculators['gdt']._cutoffs
    }
  


  def compute_fitness(self, metrics: Dict[str, float]) -> float:
    return 2.0 * metrics['gdt'] / metrics['rmsd'] + \
      (1.0 / metrics['energyScore'])



  def _sigmoid(self, x: float) -> float:
    return 1.0 / (1.0 + math.exp(x))
