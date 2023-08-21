from .FitnessFunction import FitnessFunction
from typing import Dict
from ..Metrics import EnergyScore





class InverseEnergyScore(FitnessFunction):
  
  @classmethod
  def get_name(cls) -> str:
    return 'Fitness_InverseEnergyScore'
  


  @classmethod
  def upper_bound(cls) -> float:
    return 0.95
  


  def __init__(self) -> None:
    super().__init__({ 'energy_score': EnergyScore() })
  


  def compute_fitness(self, metrics: Dict[str, float]) -> float:
    return 1.0 / metrics['energy_score']
