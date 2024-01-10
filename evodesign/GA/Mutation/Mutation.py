from abc import ABC, abstractmethod, abstractclassmethod
from evodesign.Individual import Individual
from evodesign.Population import Population
import evodesign.Choice as Choice





class Mutation(ABC):

  @abstractclassmethod
  def name(cls) -> str:
    raise NotImplementedError



  def __init__(self, probability: float = 1.0) -> None:
    super().__init__()
    self._activation_weights = ( probability, 1.0 - probability )
    self._probability = probability



  def params_json(self) -> dict:
    return {
      'probability': self._probability
    }


  
  @abstractmethod
  def mutate_sequence(self, sequence: str) -> str:
    raise NotImplementedError



  def __call__(self, children: Population) -> None:
    for i, child in enumerate(children):
      if Choice.flip_coin(self._activation_weights):
        children[i] = Individual(self.mutate_sequence(child.sequence))
