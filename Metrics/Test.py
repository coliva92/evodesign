from .Metric import Metric
from typing import List
from Bio.PDB.Atom import Atom
import random





class Test(Metric):
  
  def __call__(self, 
               modelBackbone: List[Atom], 
               referenceBackbone: List[Atom]
               ) -> float:
    return random.random()
