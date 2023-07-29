from .Recombination import Recombination
from typing import List
import random





class SinglePointCrossover(Recombination):
  """
  Operación de recombinación donde los hijos se obtienen al combinar el 
  segmento inicial de un padre con el segmento final del otro padre. El punto 
  intermedio donde se separan los segmentos para cada padre se elige de manera 
  aleatoria.
  """

  @classmethod
  def get_name(cls) -> str:
    return 'GA_Recombination_SinglePointCrossover'
  
  
  
  def __init__(self, probability: float = 1.0) -> None:
    super().__init__(probability)



  def create_children(self, 
                      mother: str,
                      father: str) -> List[str]:
    """
    Combina las secuencias especificadas por `mother` y `father` para generar 
    nuevas secuencias.
    """
    # suponemos que ambos padres son de la misma longitud
    n = len(mother)
    i = random.randint(0, n - 1)
    sister = mother[0:i] + father[i:]
    brother = father[0:i] + mother[i:]
    return [ sister, brother ]
