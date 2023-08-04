from .Replacement import Replacement
from evodesign.Population import Population





class Generational(Replacement):
  
  @classmethod
  def get_name(cls) -> str:
    return 'GA_Replacement_Generational'
  


  def __call__(self, 
               population: Population,
               children: Population
               ) -> Population:
    children.iterationId = population.iterationId + 1
    return children
