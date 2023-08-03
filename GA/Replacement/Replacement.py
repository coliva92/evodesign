from abc import ABC, abstractmethod
from typing import List
from evodesign import Individual





class Replacement(ABC):

  @classmethod
  @abstractmethod
  def get_name(cls) -> str:
    raise NotImplementedError



  @abstractmethod
  def __call__(self, 
               population: List[Individual],
               children: List[Individual]
               ) -> List[Individual]:
    raise NotImplementedError
