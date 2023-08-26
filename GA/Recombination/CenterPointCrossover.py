from .Recombination import Recombination
from typing import Tuple





class CenterPointCrossover(Recombination):
  
  @classmethod
  def name(cls) -> str:
    return 'GA_Recombination_CenterPointCrossover'



  def offspring_sequences(self, 
                          mother: str,
                          father: str
                          ) -> Tuple[str]:
    # suponemos que ambos padres son de la misma longitud
    n = len(mother)
    i = n / 2 if n % 2 == 0 else (n - 1) / 2
    sister = mother[0:i] + father[i:]
    brother = father[0:i] + mother[i:]
    return ( sister, brother )
