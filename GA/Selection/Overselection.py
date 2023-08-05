from .Selection import Selection
from typing import List
from evodesign import Individual
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
  


  def __call__(self, individuals: List[Individual]) -> List[Individual]:
    if Choice.flip_coin(self._weights):
      return random.sample(individuals[-self._top_size:], self._selection_size)
    return random.sample(individuals[0:-self._top_size], self._selection_size)
