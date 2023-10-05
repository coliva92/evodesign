from .Replacement import Replacement
from evodesign.Population import Population
from evodesign.Individual import Individual
from typing import List





class GenerationalElitism(Replacement):
  
  @classmethod
  def name(cls) -> str:
    return 'GA_Replacement_GenerationalElitism'
  


  def __init__(self, elitismSize: int = 1) -> None:
    super().__init__()
    self._elitismSize = elitismSize



  def elitism_size(self) -> int:
    return self._elitismSize
  

  
  def _merge(self, 
             a: List[Individual], 
             b: List[Individual]
             ) -> List[Individual]:
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
        break
      merged.append(b[j])
      j += 1
    return merged
  


  def __call__(self, 
               population: Population,
               children: Population
               ) -> Population:
    children.iteration_id = population.iteration_id + 1
    if self._elitismSize == 1 and children[-1] <= population[-1]:
      children.individuals = children[1:] + population[-1]
    if self._elitismSize > 1:
      top = self._merge(population[-self._elitismSize:],
                        children[-self._elitismSize:])
      children.individuals = children[:-self._elitismSize] + top
    return children
  