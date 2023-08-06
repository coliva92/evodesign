from .Recombination import Recombination
from typing import Tuple
import random





class TwoPointsCrossover(Recombination):

  @classmethod
  def get_name(cls) -> str:
    return 'GA_Recombination_TwoPointsCrossover'



  def create_offspring_sequences(self, 
                                 mother: str,
                                 father: str
                                 ) -> Tuple[str]:
    # suponemos que ambos padres son de la misma longitud
    n = len(mother)
    i = random.randint(0, n - 1) 
    j = random.randint(i, n - 1)
    sister = mother[0:i] + father[i:j] + mother[j:]
    brother = father[0:i] + mother[i:j] + father[j:]
    return ( sister, brother )
