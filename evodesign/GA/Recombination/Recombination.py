from abc import ABC, abstractmethod, abstractclassmethod
from typing import List
from ... import Individual, Population
from evodesign import Random as Choice 





class Recombination(ABC):

  @abstractclassmethod
  def name(cls) -> str:
    raise NotImplementedError



  def __init__(self, probability: float = 1.0) -> None:
    super().__init__()
    self._activation_weights = ( probability, 1.0 - probability )
    self._probability = probability
  


  def params_as_dict(self) -> dict:
    return {
      'probability': self._probability
    }



  def __call__(self, parents: Population) -> Population:
    # suponemos que la recombinación siempre se lleva a cabo entre pares
    if len(parents) % 2 != 0:
      parents = parents[:-1]
    
    n = len(parents[0])
    
    children = []
    for mother, father in zip(parents[0::2], parents[1::2]):
      if Choice.coin_toss(self._activation_weights):
        children += [
          Individual(sequence) 
          for sequence in self.offspring_sequences(mother.sequence, 
                                                   father.sequence)
        ]
        continue

      # suponemos que la recombinación siempre produce dos hijos
      children += [ Individual.new_random(n), Individual.new_random(n) ]
    return Population(children)



  @abstractmethod
  def offspring_sequences(self, 
                          mother: str,
                          father: str
                          ) -> List[str]:
    raise NotImplementedError
  