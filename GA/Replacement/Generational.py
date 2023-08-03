from .Replacement import Replacement
from typing import List
from evodesign import Individual





class Generational(Replacement):
  
  @classmethod
  def get_name(cls) -> str:
    return 'GA_Replacement_Generational'
  


  def __call__(self, 
               population: List[Individual],
               children: List[Individual]
               ) -> List[Individual]:
    return sorted(children)
