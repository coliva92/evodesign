from .Recombination import Recombination
from typing import Tuple
import random





class UniformCrossover(Recombination):
  
  _options = [ 0, 1 ]



  @classmethod
  def get_name(cls) -> str:
    return 'GA_Recombination_UniformCrossover'
  
  

  def __init__(self, 
               probability: float = 1.0, 
               bias: float = 0.5) -> None:
    super().__init__(probability)
    self._weights = [ bias, 1.0 - bias ]
    self._bias = bias
  


  def as_dict(self) -> dict:
    params = super().as_dict()
    params['bias'] = self._bias
    return params



  def create_offspring_sequences(self, 
                                 mother: str,
                                 father: str
                                 ) -> Tuple[str]:
    # suponemos que ambos padres son de la misma longitud y que vienen 
    # ordenados por aptitud de manera ascendente
    selections = random.choices(UniformCrossover._options, 
                                self._weights, 
                                k=len(mother))
    temp = ( mother, father )
    sister = ''.join([ temp[parent][i] for i, parent in enumerate(selections) ])
    temp = ( father, mother )
    brother = ''.join(
      [ temp[parent][i] for i, parent in enumerate(selections) ])
    return ( sister, brother )
