from .Mutation import Mutation
import random
from evodesign.Sequence import AMINOACIDS





class MultipleSwitches(Mutation):
  """
  Operación de mutación donde se elige aleatoriamente un residuo en la 
  secuencia y se cambia por otro aminoácido aleatorio.
  """

  def __init__(self, probability: float = 1.0, numSwitches: int = 1) -> None:
    """
    Constructor.
    - `probability`: la probabilidad con la que se aplicará la mutación sobre 
      un individuo.
    """
    super().__init__(probability)
    self._num_switches = numSwitches
  


  @classmethod
  def get_name(cls) -> str:
    return 'GA_Mutation_MultipleSwitches'



  def mutate_sequence(self, sequence: str) -> str:
    """
    Muta la secuencia especificada por `sequence` y retorna la nueva secuencia.
    """
    positions = random.sample(range(len(sequence)), self._num_switches)
    # en Python, las cadenas son inmutables; hay que convertirlas a una lista
    temp = list(sequence)
    for i in positions:
      letter = temp[i]
      while temp[i] == letter:
        letter = random.choices(AMINOACIDS, k=1)[0]
      temp[i] = letter
    return ''.join(temp)
