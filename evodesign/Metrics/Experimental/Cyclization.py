from ..Metric import Metric
from typing import List
from Bio.PDB.Atom import Atom





class Cyclization(Metric):

  @classmethod
  def name(cls) -> str:
    return 'Cyclization'



  def __call__(self, 
               sequence: str,
               model: List[Atom], 
               reference: List[Atom]
               ) -> float:
    # distance between the N atom of the first residue and the C 
    # atom of the last
    return model[0] - model[-2]
