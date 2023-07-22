from .Recombination import Recombination
from typing import List
import random





class UniformCrossover(Recombination):
  """
  Operación de recombinación donde los hijos se obtienen al combinar un 
  segmento intermedio de un padre con los segmentos extremos del otro padre. 
  Los dos puntos intermedios que separan las secuencias de cada padre en tres 
  segmentos se eligen de manera aleatoria. 
  """
  
  _options = [ 0, 1 ]

  def __init__(self, 
               probability: float = 1.0, 
               bias: float = 0.5) -> None:
    super().__init__(probability)
    self._weights = [ bias, 1.0 - bias ]
    self._bias = bias



  @classmethod
  def get_name(cls) -> str:
    return 'GA_Recombination_UniformCrossover'
  


  def get_params_memento(self) -> dict:
    params = super().get_params_memento()
    params['bias'] = self._bias
    return params



  def create_children(self, 
                      mother: str,
                      father: str) -> List[str]:
    """
    Combina las secuencias especificadas por `mother` y `father` para generar 
    nuevas secuencias.
    """
    # suponemos que ambos padres son de la misma longitud y que vienen 
    # ordenados por aptitud de manera ascendente
    n = len(mother)
    selections = random.choices(self._options, self._weights, k=n)
    temp = [ mother, father ]
    sister = ''.join([ temp[parent][i] for i, parent in enumerate(selections) ])
    temp[0], temp[1] = father, mother
    brother = ''.join([ temp[parent][i] for i, parent in enumerate(selections) ])
    return [ sister, brother ]
