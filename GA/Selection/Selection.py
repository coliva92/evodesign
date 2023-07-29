from abc import ABC, abstractmethod
from typing import List
from evodesign import Individual





class Selection(ABC):
  """
  La representación de la operación de selección de padres que se aplica para 
  seleccionar un subconjunto de individuos (i.e., de secuencias) de una 
  población particular.
  """

  @classmethod
  @abstractmethod
  def get_name(cls) -> str:
    raise NotImplementedError
  
  
  
  def __init__(self, selectionSize: int) -> None:
    """
    Constructor.
    - `selectionSize`: el número de individuos a elegir como padres.
    """
    super().__init__()
    self._selection_size = selectionSize



  @abstractmethod
  def select_parents(self, population: List[Individual]) -> List[Individual]:
    """
    Selecciona algunos individuos del arreglo especificado por `population` 
    para posteriormente recombinarlos y generar la población de la siguiente 
    generación.
    """
    raise NotImplementedError



  def get_params_memento(self) -> dict:
    return {
      'selectionSize': self._selection_size
    }



  def apply(self, population: List[Individual]) -> List[Individual]:
    return self.select_parents(population)
