from .Replacement import Replacement
from typing import List
from evodesign.Population import Individual





class WorseOut(Replacement):
  """
  Operación de reemplazo donde los peores individuos (i.e., aquellos con el 
  valor de aptitud menos favorable) son excluidos de la población de la 
  siguiente generación.
  """
  
  def __init__(self) -> None:
    super().__init__()
  


  @classmethod
  def get_name(cls) -> str:
    return 'GA_Replacement_WorseOut'



  def apply(self, 
            population: List[Individual],
            children: List[Individual]) -> List[Individual]:
    """
    Reemplaza parte de los individuos especificados por `population` con 
    aquellos especificados por `children`.
    """
    # suponemos que ambos arreglos ya tienen su valores de aptitud calculados
    next_population = sorted(population + children)
    # suponemos que estamos maximizando la aptitud
    return next_population[len(children):]
