from .Predictor import Predictor
from typing import List
from Bio.PDB.Atom import Atom
from ..Exceptions import HttpInternalServerError
import evodesign.Choice as Choice





class RemoteApiTest(Predictor):

  @classmethod
  def name(cls) -> str:
    return 'Predictor_RemoteApiTest'
  

  
  def predict_structure(self, 
                        sequence: str, 
                        pdbFilename: str
                        ) -> None:
    if Choice.flip_coin(( 0.05, 0.95 )):
      raise HttpInternalServerError
  


  def __call__(self,
               sequence: str, 
               pdbFilename: str
               ) -> List[Atom]:
    self.predict_structure(sequence, pdbFilename)
    return []
