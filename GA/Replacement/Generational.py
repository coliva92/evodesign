from .Replacement import Replacement
from evodesign.Population import Population





class Generational(Replacement):
  
  @classmethod
  def name(cls) -> str:
    return 'GA_Replacement_Generational'



  def __call__(self, 
               population: Population,
               children: Population
               ) -> Population:
    children.iteration_id = population.iteration_id + 1
    return children
  