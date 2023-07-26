from .Replacement import Replacement
from typing import List
from evodesign.Population import Individual





class Generational(Replacement):
  
  def __init__(self) -> None:
    super().__init__()
  


  @classmethod
  def get_name(cls) -> str:
    return 'GA_Replacement_Generational'



  def apply(self, 
            population: List[Individual],
            children: List[Individual]) -> List[Individual]:
    return sorted(children)
