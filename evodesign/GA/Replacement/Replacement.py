from abc import ABC, abstractmethod, abstractclassmethod
from ... import Population





class Replacement(ABC):

  @abstractclassmethod
  def name(cls) -> str:
    raise NotImplementedError



  @abstractmethod
  def __call__(self, 
               population: Population,
               children: Population
               ) -> Population:
    raise NotImplementedError
