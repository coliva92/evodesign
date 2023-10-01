from .Replacement import Replacement
from typing import List
from evodesign.Population import Population
from Individual import Individual





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
    return children



  def _merge(a: List[Individual], b: List[Individual]) -> List[Individual]:
    merged, i, j = [], 0, 0
    while True:
      if i == len(a):
        merged += b[j:]
        break
      if j == len(b):
        merged += a[i:]
        break
      if a[i] <= b[j]:
        merged.append(a[i])
        i += 1
      else:
        merged.append(b[j])
        j += 1
    return merged
  