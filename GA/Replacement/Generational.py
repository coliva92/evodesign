from .Replacement import Replacement
from typing import List
from evodesign.Population import Population





class Generational(Replacement):
  
  @classmethod
  def name(cls) -> str:
    return 'GA_Replacement_Generational'
  


  def __init__(self, elitismSize: int) -> None:
    super().__init__()
    self._elitismSize = elitismSize
  


  def __call__(self, 
               population: Population,
               children: Population
               ) -> Population:
    children.iteration_id = population.iteration_id + 1
    children.individuals = children.individuals[:-self._elitismSize] + \
      population.individuals[-self._elitismSize:]
    return children
