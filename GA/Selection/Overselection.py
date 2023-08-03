from .Selection import Selection
from typing import List
from evodesign import Individual
import random





class Overselection(Selection):

  _options = [ True, False ]



  @classmethod
  def get_name(cls) -> str:
    return 'GA_Selection_Overselection'
  


  def __init__(self, 
               selectionSize: int,
               topSize: int,
               topProbability: float = 0.8) -> None:
    super().__init__(selectionSize)
    self._top_size = topSize
    self._weights = [ topProbability, 1.0 - topProbability ]
    self._top_probability = topProbability
  


  def as_dict(self) -> dict:
    params = super().as_dict()
    params['topSize'] = self._top_size
    params['topProbability'] = self._top_probability
    return params
  


  def __call__(self, population: List[Individual]) -> List[Individual]:
    if random.choices(Overselection._options, self._weights, k=1)[0]:
      return random.sample(population[-self._top_size:], self._selection_size)
    return random.sample(population[0:-self._top_size], self._selection_size)
