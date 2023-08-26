from .FitnessFunction import FitnessFunction
from typing import Dict
from ..Metrics import Rmsd





class NegativeRmsd(FitnessFunction):
  
  @classmethod
  def name(cls) -> str:
    return 'Fitness_NegativeRMSD'
  


  @classmethod
  def upper_bound(cls) -> float:
    return -2.0
  

  
  def __init__(self) -> None:
    super().__init__({ 'rmsd': Rmsd() })
  


  def compute_fitness(self, metrics: Dict[str, float]) -> float:
    return -metrics['rmsd']
