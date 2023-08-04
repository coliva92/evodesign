from abc import ABC, abstractmethod
from evodesign.Statistics import Statistics
from evodesign.Population import Population





class Terminator(ABC):
  
  @abstractmethod
  def __call__(self,
               population: Population,
               stats: Statistics
               ) -> bool:
    raise NotImplementedError
  