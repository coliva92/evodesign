from .Recombination import Recombination
from typing import List
import random





class UniformCrossover(Recombination):

  @classmethod
  def name(cls) -> str:
    return 'GA_Recombination_UniformCrossover'
  
  

  def __init__(self, 
               probability: float = 1.0, 
               bias: float = 0.5) -> None:
    super().__init__(2, probability)
    self._weights = ( bias, 1.0 - bias )
    self._bias = bias
  


  def params_json(self) -> dict:
    params = super().params_json()
    params['bias'] = self._bias
    return params



  def offspring_sequences(self, 
                          mother: str,
                          father: str
                          ) -> List[str]:
    # suponemos que ambos padres son de la misma longitud y que vienen 
    # ordenados por aptitud de manera ascendente
    selections = random.choices(( 0, 1 ), self._weights, k=len(mother))
    parents = ( mother, father )
    sister = ''.join([ parents[p][i] for i, p in enumerate(selections) ])
    parents = ( father, mother )
    brother = ''.join([ parents[p][i] for i, p in enumerate(selections) ])
    return [ sister, brother ]
