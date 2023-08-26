from .FitnessFunction import FitnessFunction
from typing import Dict
from ..Metrics import Test as TestMetric





class Test(FitnessFunction):
  
  @classmethod
  def name(cls) -> str:
    return 'Fitness_Test'
  


  @classmethod
  def upper_bound(cls) -> float:
    return 1.0
  


  def __init__(self) -> None:
    super().__init__({ 'test': TestMetric() })
  


  def compute_fitness(self, metrics: Dict[str, float]) -> float:
    return metrics['test']
