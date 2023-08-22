from .FitnessFunction import FitnessFunction
from typing import Dict
from ..Metrics import Rmsd, Gdt as GdtMetric





class Gdt(FitnessFunction):
  
  @classmethod
  def name(cls) -> str:
    return 'Fitness_GDT'
  


  @classmethod
  def upper_bound(cls) -> float:
    return 0.95
  

  
  def __init__(self) -> None:
    super().__init__({ 'rmsd': Rmsd(), 'gdt_ts': GdtMetric() })
  


  def compute_fitness(self, metrics: Dict[str, float]) -> float:
    return metrics['gdt_ts']
