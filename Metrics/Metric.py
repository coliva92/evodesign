from abc import ABC, abstractmethod
from typing import List
from Bio.PDB.Atom import Atom





class Metric(ABC):
  """
  Una métrica de la calidad de una estructura de una proteína con respecto a 
  una estructura de referencia. 
  """
  
  def __init__(self) -> None:
    super().__init__()
  


  @abstractmethod
  def compute(self, 
              model_backbone: List[Atom], 
              reference_backbone: List[Atom]) -> float:
    """
    Calcula la métrica de calidad del esqueleto especificado por `backbone`.
    """
    pass
