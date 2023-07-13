from .Replacement import Replacement
from typing import List
from evodesign.Population import Individual
import random





class TotalReplacement(Replacement):
  """
  Operación de reemplazo donde los peores individuos (i.e., aquellos con el 
  valor de aptitud menos favorable) son excluidos de la población de la 
  siguiente generación.
  """
  
  def __init__(self) -> None:
    super().__init__()
  


  @classmethod
  def get_name(cls) -> str:
    return 'GA_Replacement_TotalReplacement'



  def apply(self, 
            population: List[Individual],
            children: List[Individual]) -> List[Individual]:
    """
    Reemplaza parte de los individuos especificados por `population` con 
    aquellos especificados por `children`.
    """
    # suponemos que el arreglo de hijos es de mayor longitud que el arreglo de
    # población
    children = sorted(children)
    m, n = len(population), len(children)
    if population[-1].fitness > children[-1].fitness:
      i = random.randint(0, n - 1)
      del children[i]
      children[-1] = population[-1]
    return children[n - m:]
