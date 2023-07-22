from .Selection import Selection
from typing import List
from evodesign.Population import Individual
import random





class Tournament(Selection):
  """
  Operación de selección donde primero el padre es el mejor individuo de un 
  grupo elegido aleatoriamente de la población. Esta operación es útil cuando 
  se tiene una población grande y no necesariamente se conoce la información de 
  la población entera.
  """
  
  def __init__(self,
               selectionSize: int,
               tournamentSize: int) -> None:
    """
    Constructor.
    - `selectionSize`: el número de individuos a elegir como padres.
    - `tournamentSize`: el número de individuos a elegir para realizar el 
      torneo.
    """
    super().__init__(selectionSize)
    self._tournament_size = tournamentSize
  


  @classmethod
  def get_name(cls) -> str:
    return 'GA_Selection_Tournament'
  


  def get_params_memento(self) -> dict:
    params = super().get_params_memento()
    params['tournamentSize'] = self._tournament_size
    return params



  def select_parents(self, population: List[Individual]) -> List[Individual]:
    """
    Selecciona algunos individuos del arreglo especificado por `population` 
    para posteriormente recombinarlos y generar la población de la siguiente 
    generación.
    """
    selected_parents = []
    for _ in range(self._selection_size):
      candidates = random.sample(population, self._tournament_size)
      candidates = sorted(candidates)
      selected_parents.append(candidates[-1])
    return selected_parents
    