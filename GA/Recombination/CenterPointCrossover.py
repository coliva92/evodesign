from .Recombination import Recombination
from typing import List





class CenterPointCrossover(Recombination):
  
  @classmethod
  def name(cls) -> str:
    return 'GA_Recombination_CenterPointCrossover'
  


  def __init__(self, probability: float = 1.0) -> None:
    super().__init__(2, probability)



  def offspring_sequences(self, 
                          mother: str,
                          father: str
                          ) -> List[str]:
    # suponemos que ambos padres son de la misma longitud
    n = len(mother)
    i = n / 2 if n % 2 == 0 else (n - 1) / 2
    sister = mother[0:i] + father[i:]
    brother = father[0:i] + mother[i:]
    return [ sister, brother ]
