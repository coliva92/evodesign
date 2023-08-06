from abc import ABC, abstractmethod
from typing import List
from evodesign.Population import Population
from evodesign.Individual import Individual





class Selection(ABC):

  @classmethod
  @abstractmethod
  def get_name(cls) -> str:
    raise NotImplementedError
  
  
  
  def __init__(self, selectionSize: int) -> None:
    super().__init__()
    self._selection_size = selectionSize



  def as_json(self) -> dict:
    return {
      'selectionSize': self._selection_size
    }



  @abstractmethod
  def select_parents(self, population: Population) -> List[Individual]:
    raise NotImplementedError
  


  def __call__(self, population: Population) -> Population:
    return Population(self.select_parents(population))  
