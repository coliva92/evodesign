from .Metric import Metric
from Bio.PDB.Atom import Atom
from typing import List
import random





class Testing(Metric):
  
  def __init__(self) -> None:
    super().__init__()
  


  def compute(self, 
              model_backbone: List[Atom],
              reference_backbone: List[Atom]) -> float:
    """
    Superpone el esqueleto especificado por `backbone` sobre el esqueleto de 
    referencia y retorna el RMSD resultante.
    """
    return random.random()
  