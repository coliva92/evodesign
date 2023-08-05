from .Selection import Selection
from typing import List
from evodesign.Population import Population
from evodesign.Individual import Individual
import random





class Uniform(Selection):

  @classmethod
  def get_name(cls) -> str:
    return 'GA_Selection_Uniform'



  def __call__(self, population: Population) -> List[Individual]:
    return random.sample(population.individual, self._selection_size)
