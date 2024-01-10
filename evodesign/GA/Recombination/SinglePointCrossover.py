from .Recombination import Recombination
from typing import List
import random





class SinglePointCrossover(Recombination):

  @classmethod
  def name(cls) -> str:
    return 'GA_Recombination_SinglePointCrossover'
  


  def __init__(self, probability: float = 1) -> None:
    super().__init__(2, probability)



  def offspring_sequences(self, 
                          mother: str,
                          father: str
                          ) -> List[str]:
    # suponemos que ambos padres son de la misma longitud
    i = random.randint(0, len(mother) - 1)
    sister = mother[0:i] + father[i:]
    brother = father[0:i] + mother[i:]
    return [ sister, brother ]
