from . import Selection
from typing import List
from ... import Population
from ... import Individual
import random





class Uniform(Selection):

  @classmethod
  def name(cls) -> str:
    return 'GA_Selection_Uniform'



  def select_parents(self, population: Population) -> List[Individual]:
    return random.sample(population.individuals, self._selection_size)
