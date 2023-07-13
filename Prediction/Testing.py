from .Predictor import Predictor
from typing import List
from Bio.PDB.Atom import Atom
import random





class Testing(Predictor):
  
  def __init__(self) -> None:
    super().__init__()
  


  @classmethod
  def get_name(cls) -> str:
    return 'Predictor_Testing'
  


  def get_predicted_backbone(self, 
                             id: str,
                             sequence: str, 
                             pdbFilename: str) -> List[Atom]:
    return []
  
  
    
  def predict_structure(self, 
                        sequence: str, 
                        pdbFilename: str) -> None:
    """
    Predice la estructura de la secuencia de aminiácidos especificada por 
    `sequence` y escribe en resultado en un archivo PDB cuyo nombre está 
    especificado por `pdbFilename`.
    """
    if random.randint(0, 1000) < 10:
      raise RuntimeError('simulating a failed call to a remote API')
