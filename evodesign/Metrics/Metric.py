from abc import ABC, abstractmethod
from typing import List
from Bio.PDB.Atom import Atom





class Metric(ABC):
  
  @abstractmethod
  def __call__(self, 
               sequence: str,
               model: List[Atom], 
               reference: List[Atom]
               ) -> float:
    raise NotImplementedError
