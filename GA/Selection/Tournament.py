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
  


  def as_dict(self) -> dict:
    params = super().as_dict()
    params['tournamentSize'] = self._tournament_size
    return params



  def __call__(self, population: List[Individual]) -> List[Individual]:
    selected_parents = []
    for i in range(self._selection_size):
      candidates = random.sample(population, self._tournament_size)
      winner = sorted(candidates)[-1]
      # garantizamos que dos padres consecutivos siempre sean diferentes
      if i % 2 != 0:
        while selected_parents[i - 1].sequence == winner.sequence:
          candidates = random.sample(population, self._tournament_size)
          winner = sorted(candidates)[-1]
      selected_parents.append(winner)
    return selected_parents
    