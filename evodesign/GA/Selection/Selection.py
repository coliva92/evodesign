from abc import ABC, abstractmethod, abstractclassmethod
from typing import List
from ... import Population
from ... import Individual





class Selection(ABC):

  @abstractclassmethod
  def name(cls) -> str:
    raise NotImplementedError
  
  
  
  def __init__(self, selectionSize: int) -> None:
    super().__init__()
    self._selection_size = selectionSize



  def params_as_dict(self) -> dict:
    return {}



  @abstractmethod
  def select_parents(self, population: Population) -> List[Individual]:
    raise NotImplementedError
  


  def __call__(self, population: Population) -> Population:
    return Population(self.select_parents(population))  
