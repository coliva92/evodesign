from abc import ABC, abstractmethod
from typing import List
from evodesign.Population import Population
from evodesign.Individual import Individual





class Selection(ABC):

  @classmethod
  @abstractmethod
  def name(cls) -> str:
    raise NotImplementedError
  
  
  
  def __init__(self, selectionSize: int) -> None:
    super().__init__()
    # vamos a suponer que en toda operacion de recombinacion se preserva 
    # únicamente el mejor hijo y se descarta el otro; entonces, el valor de
    # `selectionSize` debe duplicarse, de lo contrario solo se seleccionarían
    # la mitad de los individuos deseados
    self._selection_size = selectionSize
    self._num_individual_selections = 2 * selectionSize



  def params_json(self) -> dict:
    return {
      'selectionSize': self._selection_size
    }



  @abstractmethod
  def select_parents(self, population: Population) -> List[Individual]:
    raise NotImplementedError
  


  def __call__(self, population: Population) -> Population:
    return Population(self.select_parents(population))  
