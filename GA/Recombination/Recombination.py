from abc import ABC, abstractmethod, abstractclassmethod
from typing import List
from evodesign import Individual, Population
import evodesign.Choice as Choice





class Recombination(ABC):

  @abstractclassmethod
  def name(cls) -> str:
    raise NotImplementedError



  def __init__(self, 
               numOffspringsPerPair: int,
               probability: float = 1.0
               ) -> None:
    super().__init__()
    self._activation_weights = ( probability, 1.0 - probability )
    self._probability = probability
    self._num_offsprings = numOffspringsPerPair
  


  def params_json(self) -> dict:
    return {
      'probability': self._probability
    }



  def __call__(self, parents: Population) -> Population:
    # suponemos que la recombinación siempre se lleva a cabo entre pares
    if len(parents) % 2 != 0:
      parents = parents[:-1]
    
    children = []
    for mother, father in zip(parents[0::2], parents[1::2]):
      if Choice.flip_coin(self._activation_weights):
        children += self.offspring_sequences(mother.sequence,
                                             father.sequence)
        continue
      children += [ 
        Individual.new_random(len(mother)) \
        for _ in range(self._num_offsprings) 
      ]
    return Population(children)



  @abstractmethod
  def offspring_sequences(self, 
                          mother: str,
                          father: str
                          ) -> List[str]:
    raise NotImplementedError
  