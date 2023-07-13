from abc import ABC, abstractmethod
import random
from typing import List
from evodesign.Population import Individual





class Mutation(ABC):
  """
  La representación para la operación de mutación que se aplica sobre los 
  individuos (i.e., sobre las secuencias de aminoácidos) de la población.
  """
  
  _options = [ True, False ]



  def __init__(self, probability: float = 1.0) -> None:
    """
    Constructor.
    - `probability`: la probabilidad con la que se aplicará la mutación sobre 
      un individuo.
    """
    super().__init__()
    self._weights = [ probability, 1.0 - probability ]
    self._probability = probability
  


  @classmethod
  @abstractmethod
  def get_name(cls) -> str:
    pass



  def get_params_memento(self) -> dict:
    return {
      'probability': self._probability
    }


  
  @abstractmethod
  def mutate_sequence(self, sequence: str) -> str:
    """
    Muta la secuencia especificada por `sequence` y retorna la nueva secuencia.
    """
    pass



  def apply(self, children: List[Individual]) -> None:
    """
    Aplica la mutación a los individuos especificados por `children` con la 
    probabilidad especificada en el constructor.
    """
    for child in children:
      if random.choices(Mutation._options, self._weights, k=1)[0]:
        child.sequence = self.mutate_sequence(child.sequence)
