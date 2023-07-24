from .Selection import Selection
from typing import List
from evodesign.Population import Individual
import random





class Uniform(Selection):
  """
  Operación de selección donde los padres se eligen aleatoriamente de toda la 
  población con probabilidad uniforme.
  """
  
  def __init__(self, selectionSize: int) -> None:
    """
    Constructor.
    - `selectionSize`: el número de individuos a elegir como padres.
    """
    super().__init__(selectionSize)
  


  @classmethod
  def get_name(cls) -> str:
    return 'GA_Selection_Uniform'



  def select_parents(self, population: List[Individual]) -> List[Individual]:
    """
    Selecciona algunos individuos del arreglo especificado por `population` 
    para posteriormente recombinarlos y generar la población de la siguiente 
    generación.
    """
    return random.sample(population, self._selection_size)
