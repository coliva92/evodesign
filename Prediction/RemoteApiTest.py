from .Predictor import Predictor
from typing import List
from Bio.PDB.Atom import Atom
import random
from ..Exceptions import RemoteApiRequestsExceeded





class RemoteApiTest(Predictor):

  _options = [ True, False ]
  _weights = [ 0.05, 0.95 ]



  @classmethod
  def get_name(cls) -> str:
    return 'Predictor_RemoteApiTest'
  

  
  def predict_structure(self, 
                        sequence: str, 
                        pdbFilename: str
                        ) -> None:
    if random.choices(RemoteApiTest._options, k=1)[0]:
      raise RemoteApiRequestsExceeded
  


  def __call__(self,
               sequence: str, 
               pdbFilename: str
               ) -> List[Atom]:
    self.predict_structure(sequence, pdbFilename)
    return []
