from ..Metric import Metric
from typing import List, Optional
from Bio.PDB.Atom import Atom





class Cyclization(Metric):

  @classmethod
  def column_name(cls) -> str:
    return 'Cyclization'



  def __call__(self, 
               model: List[Atom], 
               reference: List[Atom],
               sequence: Optional[str] = None
               ) -> float:
    # distance between the N atom of the first residue and the C 
    # atom of the last
    return model[0] - model[-2]
