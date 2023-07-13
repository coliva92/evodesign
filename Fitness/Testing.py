from .FitnessFunction import FitnessFunction
from ..Metrics import Testing as TestMetric
from typing import Dict





class Testing(FitnessFunction):
  
  def __init__(self) -> None:
    """
    Constructor.
    """
    super().__init__({ 'value': TestMetric() })
  


  @classmethod
  def get_name(cls) -> str:
    return 'Fitness_Testing'
  


  @classmethod
  def upper_bound(cls) -> float:
    return 0.0
  


  def compute_fitness(self, metrics: Dict[str, float]) -> float:
    """
    Calcula la aptitud de la secuencia especificada por `sequence`, utilizando
    los valores especificados por `metrics`.
    """
    return 0.0 - metrics['value']
