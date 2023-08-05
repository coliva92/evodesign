from .Predictor import Predictor
from typing import List
from Bio.PDB.Atom import Atom
from ..Exceptions import RemoteApiRequestsExceeded
import evodesign.Choice as Choice





class RemoteApiTest(Predictor):

  @classmethod
  def get_name(cls) -> str:
    return 'Predictor_RemoteApiTest'
  

  
  def predict_structure(self, 
                        sequence: str, 
                        pdbFilename: str
                        ) -> None:
    if Choice.flip_coin(( 0.05, 0.95 )):
      raise RemoteApiRequestsExceeded
  


  def __call__(self,
               sequence: str, 
               pdbFilename: str
               ) -> List[Atom]:
    self.predict_structure(sequence, pdbFilename)
    return []
