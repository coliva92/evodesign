from abc import ABC, abstractmethod
from ... import Statistics
from ... import Population





class Terminator(ABC):
  
  @abstractmethod
  def __call__(self,
               population: Population,
               stats: Statistics
               ) -> bool:
    raise NotImplementedError
  