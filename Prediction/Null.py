from .Predictor import Predictor
from typing import List
from Bio.PDB.Atom import Atom





class Null(Predictor):
  
  @classmethod
  def name(cls) -> str:
    return 'Predictor_Null'
  


  def predict_structure(self, 
                        sequence: str, 
                        pdbFilename: str) -> None:
    pass



  def __call__(self,
               sequence: str, 
               pdbFilename: str
               ) -> List[Atom]:
    return []
