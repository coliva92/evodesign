from abc import ABC, abstractmethod
from typing import List, Tuple
from evodesign import Individual
from evodesign.Sequence import create_random_sequence
import random





class Recombination(ABC):

  @classmethod
  @abstractmethod
  def get_name(cls) -> str:
    raise NotImplementedError



  def __init__(self, probability: float = 1.0) -> None:
    super().__init__()
    self._weights = [ probability, 1.0 - probability ]
    self._probability = probability
  


  def as_json(self) -> dict:
    return {
      'probability': self._probability
    }


  
  @abstractmethod
  def create_offspring_sequences(self, 
                                 mother: str,
                                 father: str
                                 ) -> Tuple[str]:
    raise NotImplementedError



  def __call__(self, parents: List[Individual]) -> List[Individual]:
    # suponemos que la recombinación siempre se lleva a cabo entre pares y que,
    # por ende, siempre produce dos hijos
    if len(parents) % 2 != 0:
      parents = parents[:-1]
    
    def offspring(flag: bool, 
                  mother: str, 
                  father: str
                  ) -> Tuple[str]:
      return self.create_offspring_sequences(mother.sequence, father.sequence) \
        if flag else ( create_random_sequence(len(mother)), 
                       create_random_sequence(len(mother)) )
    
    flags = random.choices([ True, False ], 
                           self._weights, 
                           k=int(len(parents) / 2))
    return [
      Individual(sequence)
      for flag, mother, father in zip(flags, parents[0::2], parents[1::2])
      for sequence in offspring(flag, mother, father)
    ]
