from abc import ABC, abstractmethod
import random
from typing import List
from evodesign.Individual import Individual
import evodesign.Choice as Choice





class Mutation(ABC):

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
  def mutate_sequence(self, sequence: str) -> str:
    raise NotImplementedError



  def __call__(self, children: List[Individual]) -> None:
    for i, child in enumerate(children):
      if Choice.flip_coin(self._weights):
        children[i] = Individual(self.mutate_sequence(child.sequence))
