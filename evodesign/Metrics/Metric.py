from abc import ABC, abstractmethod, abstractclassmethod
from typing import List, Optional
from Bio.PDB.Atom import Atom





class Metric(ABC):

  @abstractclassmethod
  def name(cls) -> str:
    raise NotImplemented
  

  
  @abstractmethod
  def __call__(self, 
               model: List[Atom], 
               reference: List[Atom],
               sequence: Optional[str] = None
               ) -> float:
    raise NotImplementedError
