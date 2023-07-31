from abc import ABC, abstractmethod
import random
from typing import List
from evodesign import Individual





class Mutation(ABC):
  """
  La representación para la operación de mutación que se aplica sobre los 
  individuos (i.e., sobre las secuencias de aminoácidos) de la población.
  """
  
  _options = [ True, False ]



  @classmethod
  @abstractmethod
  def get_name(cls) -> str:
    raise NotImplementedError



  def __init__(self, probability: float = 1.0) -> None:
    """
    Constructor.
    - `probability`: la probabilidad con la que se aplicará la mutación sobre 
      un individuo.
    """
    super().__init__()
    self._weights = [ probability, 1.0 - probability ]
    self._probability = probability



  def get_params_memento(self) -> dict:
    return {
      'probability': self._probability
    }


  
  @abstractmethod
  def mutate_sequence(self, sequence: str) -> str:
    """
    Muta la secuencia especificada por `sequence` y retorna la nueva secuencia.
    """
    raise NotImplementedError



  def apply(self, children: List[Individual]) -> None:
    """
    Aplica la mutación a los individuos especificados por `children` con la 
    probabilidad especificada en el constructor.
    """
    for i, child in enumerate(children):
      if random.choices(Mutation._options, self._weights, k=1)[0]:
        children[i] = Individual(self.mutate_sequence(child.sequence))
