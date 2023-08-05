from abc import ABC, abstractmethod
from typing import List
from evodesign.Individual import Individual
from evodesign.Population import Population





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
  def __call__(self, population: Population) -> List[Individual]:
    raise NotImplementedError
