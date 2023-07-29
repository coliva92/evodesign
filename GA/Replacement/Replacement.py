from abc import ABC, abstractmethod
from typing import List
from evodesign.Population import Individual





class Replacement(ABC):
  """
  La representación para la operación de reemplazo que se aplica sobre 
  algunos individuos de la población (i.e., sobre algunas secuencias).
  """
  
  def __init__(self) -> None:
    super().__init__()
  


  @classmethod
  @abstractmethod
  def get_name(cls) -> str:
    raise NotImplementedError



  def get_params_memento(self) -> dict:
    return {}



  @abstractmethod
  def apply(self, 
            population: List[Individual],
            children: List[Individual]) -> List[Individual]:
    """
    Reemplaza parte de los individuos especificados por `population` con 
    aquellos especificados por `children`.
    """
    raise NotImplementedError
