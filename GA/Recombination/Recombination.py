from abc import ABC, abstractmethod
from typing import List
from evodesign import Individual
import random





class Recombination(ABC):
  """
  La representación para la operación de recombinación que se aplica sobre 
  algunos individuos de la población (i.e., sobre algunas secuencias).
  """

  _options = [ True, False ]



  @classmethod
  @abstractmethod
  def get_name(cls) -> str:
    raise NotImplementedError



  def __init__(self, probability: float = 1.0) -> None:
    """
    Constructor.
    - `probability`: la probabilidad con la que se desea aplicar la 
    recombinación.
    """
    super().__init__()
    self._weights = [ probability, 1.0 - probability ]
    self._probability = probability
  


  def get_params_memento(self) -> dict:
    return {
      'probability': self._probability
    }


  
  @abstractmethod
  def create_children(self, 
                      mother: str,
                      father: str) -> List[str]:
    """
    Combina las secuencias especificadas por `mother` y `father` para generar 
    nuevas secuencias.
    """
    raise NotImplementedError



  def apply(self, parents: List[Individual]) -> List[Individual]:
    """
    Combinamos los individuos especificados en `parents` para generar nuevos 
    individuos (i.e., nuevas secuencias).
    """
    # suponemos que la combinación siempre se hace en pares
    num_recombinations = len(parents)
    if num_recombinations % 2 != 0:
      num_recombinations -= 1
    children = []
    i = 0
    while (i < num_recombinations):
      if random.choices(Recombination._options, self._weights, k=1)[0]:
        sequences = self.create_children(parents[i].sequence,
                                         parents[i + 1].sequence)
        sister, brother = Individual(sequences[0]), Individual(sequences[1])
      else:
        n = len(parents[i].sequence)
        sister, brother = Individual.random(n), Individual.random(n)
      children.append(sister)
      children.append(brother)
      i += 2
    return children
