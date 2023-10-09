from .Selection import Selection
from typing import List
from evodesign.Population import Population
from evodesign.Individual import Individual
import random





class Uniform(Selection):

  @classmethod
  def name(cls) -> str:
    return 'GA_Selection_Uniform'



  def select_parents(self, population: Population) -> List[Individual]:
    return random.sample(population.individuals, self._num_individual_selections)
