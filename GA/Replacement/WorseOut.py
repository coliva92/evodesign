from .Replacement import Replacement
from typing import List
from evodesign import Individual





class WorseOut(Replacement):
  """
  Operación de reemplazo donde los peores individuos (i.e., aquellos con el 
  valor de aptitud menos favorable) son excluidos de la población de la 
  siguiente generación.
  """
  
  @classmethod
  def get_name(cls) -> str:
    return 'GA_Replacement_WorseOut'
  

  
  def __init__(self) -> None:
    super().__init__()



  def apply(self, 
            population: List[Individual],
            children: List[Individual]) -> List[Individual]:
    """
    Reemplaza parte de los individuos especificados por `population` con 
    aquellos especificados por `children`.
    """
    for child in children:
      for i, individual in enumerate(population):
        if child > individual:
          continue
        population.insert(i, child)
        break
    return population[len(children):]
