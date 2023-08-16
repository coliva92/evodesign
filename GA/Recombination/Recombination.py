from abc import ABC, abstractmethod
from typing import List, Tuple
from evodesign import Individual
from evodesign.Population import Population
import evodesign.Choice as Choice





class Recombination(ABC):

  @classmethod
  @abstractmethod
  def get_name(cls) -> str:
    raise NotImplementedError



  def __init__(self, probability: float = 1.0) -> None:
    super().__init__()
    self._activation_weights = ( probability, 1.0 - probability )
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



  def __call__(self, parents: Population) -> Population:
    # suponemos que la recombinación siempre se lleva a cabo entre pares y que,
    # por ende, siempre produce dos hijos
    if len(parents) % 2 != 0:
      parents = parents[:-1]
    
    children = []
    for mother, father in zip(parents[0::2], parents[1::2]):
      if Choice.flip_coin(self._activation_weights):
        sister, brother = self.create_offspring_sequences(mother.sequence,
                                                          father.sequence)
        children.append(Individual(sister))
        children.append(Individual(brother))
        continue
      children.append(Individual.new_random(len(mother)))
      children.append(Individual.new_random(len(mother)))
    return Population(children)
