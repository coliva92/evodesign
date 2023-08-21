from .Predictor import Predictor
from typing import List
from Bio.PDB.Atom import Atom
import os





class Null(Predictor):
  
  @classmethod
  def get_name(cls) -> str:
    return 'Predictor_Null'
  


  def predict_structure(self, 
                        sequence: str, 
                        pdbFilename: str) -> None:
    pass



  def __call__(self,
               sequence: str, 
               pdbFilename: str
               ) -> List[Atom]:
    if not os.path.isfile(pdbFilename):
      os.makedirs(os.path.dirname(os.path.abspath(pdbFilename)), exist_ok=True)
    return []
