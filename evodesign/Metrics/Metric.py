from abc import ABC, abstractmethod, abstractclassmethod
from typing import List
from Bio.PDB.Atom import Atom





class Metric(ABC):

  @abstractclassmethod
  def name(cls) -> str:
    raise NotImplemented
  

  
  @abstractmethod
  def __call__(self, 
               sequence: str,
               model: List[Atom], 
               reference: List[Atom]
               ) -> float:
    raise NotImplementedError
