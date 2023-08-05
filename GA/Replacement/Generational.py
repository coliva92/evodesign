from .Replacement import Replacement
from typing import List
from evodesign.Population import Population
from evodesign.Individual import Individual





class Generational(Replacement):
  
  @classmethod
  def get_name(cls) -> str:
    return 'GA_Replacement_Generational'
  


  def __call__(self, 
               population: Population,
               children: List[Individual]
               ) -> Population:
    children.iterationId = population.iterationId + 1
    return Population(population.iterationId + 1, children)
