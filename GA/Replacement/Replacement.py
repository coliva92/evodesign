from abc import ABC, abstractmethod
from typing import List
from evodesign.Population import Population





class Replacement(ABC):

  @classmethod
  @abstractmethod
  def get_name(cls) -> str:
    raise NotImplementedError



  @abstractmethod
  def __call__(self, 
               population: Population,
               children: Population
               ) -> Population:
    raise NotImplementedError
