from .FitnessFunction import FitnessFunction
from typing import Dict
from ..Metrics import Rmsd





class NegativeRmsd(FitnessFunction):
  """
  La función objetivo que consiste únicamente de minimizar el RMSD.
  """
  
  @classmethod
  def get_name(cls) -> str:
    return 'Fitness_NegativeRMSD'
  


  @classmethod
  def upper_bound(cls) -> float:
    return -1.95
  

  
  def __init__(self) -> None:
    super().__init__({ 'rmsd': Rmsd() })
  


  def compute_fitness(self, 
                      sequence: str, 
                      metrics: Dict[str, float]
                      ) -> float:
    """
    Calcula la aptitud de la secuencia especificada por `sequence`, utilizando
    los valores especificados por `metrics`.
    """
    return -metrics['rmsd']
