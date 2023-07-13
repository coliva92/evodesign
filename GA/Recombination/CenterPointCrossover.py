from .Recombination import Recombination
from typing import List





class CenterPointCrossover(Recombination):
  """
  Operación de recombinación donde los hijos se obtienen al combinar la mitad 
  de un padre con la mitad del otro padre.
  """
  
  def __init__(self, probability: float = 1.0) -> None:
    super().__init__(probability)
  


  @classmethod
  def get_name(cls) -> str:
    return 'GA_Recombination_CenterPointCrossover'



  def create_children(self, 
                      mother: str,
                      father: str) -> List[str]:
    """
    Combina las secuencias especificadas por `mother` y `father` para generar 
    nuevas secuencias.
    """
    # suponemos que ambos padres son de la misma longitud
    n = len(mother)
    i = int(n / 2) if n % 2 == 0 else int((n - 1) / 2)
    sister = mother[0:i] + father[i:]
    brother = father[0:i] + mother[i:]
    return [ sister, brother ]
