from .Selection import Selection
from typing import List
from evodesign import Individual
import random





class Uniform(Selection):

  @classmethod
  def get_name(cls) -> str:
    return 'GA_Selection_Uniform'



  def __call__(self, population: List[Individual]) -> List[Individual]:
    return random.sample(population, self._selection_size)
