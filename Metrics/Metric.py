from abc import ABC, abstractmethod
from typing import List
from Bio.PDB.Atom import Atom





class Metric(ABC):
  
  @abstractmethod
  def __call__(self, 
               modelBackbone: List[Atom], 
               referenceBackbone: List[Atom],
               sequence: str
               ) -> float:
    raise NotImplementedError
