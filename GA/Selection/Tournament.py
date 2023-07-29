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

  @classmethod
  def get_name(cls) -> str:
    return 'GA_Selection_Tournament'
  
  
  
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
    for i in range(self._selection_size):
      candidates = random.sample(population, self._tournament_size)
      candidates = sorted(candidates)
      winner = candidates[-1]
      # garantizamos que dos padres consecutivos siempre sean diferentes
      if i % 2 != 0:
        while selected_parents[i - 1].sequence == winner.sequence:
          candidates = random.sample(population, self._tournament_size)
          candidates = sorted(candidates)
          winner = candidates[-1]
      selected_parents.append(winner)
    return selected_parents
    