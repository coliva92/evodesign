from .Recombination import Recombination
from typing import Tuple
import random





class SinglePointCrossover(Recombination):

  @classmethod
  def name(cls) -> str:
    return 'GA_Recombination_SinglePointCrossover'



  def offspring_sequences(self, 
                          mother: str,
                          father: str
                          ) -> Tuple[str]:
    # suponemos que ambos padres son de la misma longitud
    i = random.randint(0, len(mother) - 1)
    sister = mother[0:i] + father[i:]
    brother = father[0:i] + mother[i:]
    return ( sister, brother )
