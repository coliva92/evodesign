from .Selection import Selection
from typing import List
from evodesign import Individual
import random





class Tournament(Selection):

  @classmethod
  def get_name(cls) -> str:
    return 'GA_Selection_Tournament'
  
  
  
  def __init__(self,
               selectionSize: int,
               tournamentSize: int) -> None:
    super().__init__(selectionSize)
    self._tournament_size = tournamentSize
  


  def as_json(self) -> dict:
    params = super().as_json()
    params['tournamentSize'] = self._tournament_size
    return params



  def __call__(self, individuals: List[Individual]) -> List[Individual]:
    selected_parents = []
    for i in range(self._selection_size):
      winner = sorted(random.sample(individuals, self._tournament_size))[-1]
      # garantizamos que dos padres consecutivos siempre sean diferentes
      while i % 2 != 0 and selected_parents[i - 1].sequence == winner.sequence:
        winner = sorted(random.sample(individuals, self._tournament_size))[-1]
      selected_parents.append(winner)
    return selected_parents
    