from .Replacement import Replacement
from typing import List
from evodesign import Individual





class WorseOut(Replacement):
  
  @classmethod
  def get_name(cls) -> str:
    return 'GA_Replacement_WorseOut'



  def __call__(self, 
               population: List[Individual],
               children: List[Individual]
               ) -> List[Individual]:
    for child in children:
      for i, individual in enumerate(population):
        if child > individual:
          continue
        population.insert(i, child)
        break
    return population[len(children):]
