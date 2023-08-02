from abc import ABC, abstractmethod
from typing import List
from evodesign.Statistics import Statistics
from evodesign.Individual import Individual





class Terminator(ABC):
  
  @abstractmethod
  def __call__(self, 
               iterationId: int, 
               population: List[Individual],
               stats: Statistics
               ) -> bool:
    raise NotImplementedError
  