from typing import Dict, Optional
from evodesign.Metrics import TerminalsDistance
from .FitnessFunction import FitnessFunction





class NegativeCyclization(FitnessFunction):

  @classmethod
  def name(cls) -> str:
    return 'Fitness_Cyclization'
  


  @classmethod
  def upper_bound(cls) -> float:
    return -1.4
  


  def __init__(self) -> None:
    super().__init__({
      'terminals_distance': TerminalsDistance()
    })
  


  def params_json(self) -> dict:
    return {}
  


  def compute_fitness(self, metrics: Dict[str, float], _: Optional[str] = None) -> float:
    return -metrics['terminals_distance']
  