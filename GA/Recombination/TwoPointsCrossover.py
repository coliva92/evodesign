from .Recombination import Recombination
from typing import List
import random





class TwoPointsCrossover(Recombination):

  @classmethod
  def name(cls) -> str:
    return 'GA_Recombination_TwoPointsCrossover'
  


  def __init__(self, probability: float = 1) -> None:
    super().__init__(2, probability)



  def offspring_sequences(self, 
                          mother: str,
                          father: str
                          ) -> List[str]:
    # suponemos que ambos padres son de la misma longitud
    n = len(mother)
    i = random.randint(0, n - 1) 
    j = random.randint(i, n - 1)
    sister = mother[0:i] + father[i:j] + mother[j:]
    brother = father[0:i] + mother[i:j] + father[j:]
    return [ sister, brother ]
