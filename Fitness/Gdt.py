from .FitnessFunction import FitnessFunction
from typing import Dict
from ..Metrics import Rmsd, Gdt as GdtMetric





class Gdt(FitnessFunction):
  
  def __init__(self) -> None:
    super().__init__({ 'rmsd': Rmsd(), 'gdt_ts': GdtMetric() })
  


  @classmethod
  def get_name(cls) -> str:
    return 'Fitness_GDT'
  


  @classmethod
  def upper_bound(cls) -> float:
    return 0.95
  


  def compute_fitness(self, 
                      sequence: str, 
                      metrics: Dict[str, float]
                      ) -> float:
    """
    Calcula la aptitud de la secuencia especificada por `sequence`, utilizando
    los valores especificados por `metrics`.
    """
    return metrics['gdt_ts']
