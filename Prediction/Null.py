from .Predictor import Predictor
from typing import List
from Bio.PDB.Atom import Atom





class Null(Predictor):
  
  def __init__(self) -> None:
    super().__init__()
  


  @classmethod
  def get_name(cls) -> str:
    return 'Predictor_Null'
  


  def get_predicted_backbone(self, 
                             id: str,
                             sequence: str, 
                             pdbFilename: str) -> List[Atom]:
    return []
  
  
    
  def predict_structure(self, 
                        sequence: str, 
                        pdbFilename: str) -> None:
    return
