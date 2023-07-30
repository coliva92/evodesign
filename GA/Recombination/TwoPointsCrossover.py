from .Recombination import Recombination
from typing import List
import random





class TwoPointsCrossover(Recombination):
  """
  Operación de recombinación donde los hijos se obtienen al combinar un 
  segmento intermedio de un padre con los segmentos extremos del otro padre. 
  Los dos puntos intermedios que separan las secuencias de cada padre en tres 
  segmentos se eligen de manera aleatoria. 
  """

  @classmethod
  def get_name(cls) -> str:
    return 'GA_Recombination_TwoPointsCrossover'
  
  
  
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
    j = random.randint(i, n - 1)
    sister = mother[0:i] + father[i:j] + mother[j:]
    brother = father[0:i] + mother[i:j] + father[j:]
    return [ sister, brother ]
