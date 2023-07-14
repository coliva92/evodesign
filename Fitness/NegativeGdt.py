from .FitnessFunction import FitnessFunction
from typing import Dict
from ..Metrics import Rmsd, Gdt





class NegativeGdt(FitnessFunction):
  """
  La función objetivo que consiste únicamente de minimizar el RMSD.
  """
  
  def __init__(self) -> None:
    super().__init__({ 'rmsd': Rmsd(), 'gdt_ts': Gdt() })
  


  @classmethod
  def get_name(cls) -> str:
    return 'Fitness_NegativeGdt'
  


  @classmethod
  def upper_bound(cls) -> float:
    return 0.0
  


  def compute_fitness(self, metrics: Dict[str, float]) -> float:
    """
    Calcula la aptitud de la secuencia especificada por `sequence`, utilizando
    los valores especificados por `metrics`.
    """
    return -metrics['gdt_ts']
