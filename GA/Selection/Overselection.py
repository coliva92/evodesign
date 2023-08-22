from .Selection import Selection
from typing import List
from evodesign.Individual import Individual
from evodesign.Population import Population
import random
import evodesign.Choice as Choice





class Overselection(Selection):

  @classmethod
  def get_name(cls) -> str:
    return 'GA_Selection_Overselection'
  


  def __init__(self, 
               selectionSize: int,
               topSize: int,
               topProbability: float = 0.8) -> None:
    super().__init__(selectionSize)
    self._top_size = topSize
    self._weights = ( topProbability, 1.0 - topProbability )
    self._top_probability = topProbability
  


  def as_json(self) -> dict:
    params = super().as_json()
    params['topSize'] = self._top_size
    params['topProbability'] = self._top_probability
    return params
  


  def select_parents(self, population: Population) -> List[Individual]:
    selected_parents = []
    for i in range(self._selection_size):
      pool = population[-self._top_size:] if Choice.flip_coin(self._weights) \
                                          else population[:-self._top_size]
      parent = random.choice(pool)
      # garantizamos que dos padres consecutivos siempre sean diferentes
      while i % 2 != 0 and selected_parents[i - 1].sequence == parent.sequence:
        parent = random.choice(pool)
      selected_parents.append(parent)
    return selected_parents
