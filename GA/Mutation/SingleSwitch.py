from .Mutation import Mutation
import random
from evodesign.Sequence import AMINOACIDS





class SingleSwitch(Mutation):
  """
  Operación de mutación donde se elige aleatoriamente un residuo en la 
  secuencia y se cambia por otro aminoácido aleatorio.
  """

  def __init__(self, probability: float = 1.0) -> None:
    """
    Constructor.
    - `probability`: la probabilidad con la que se aplicará la mutación sobre 
      un individuo.
    """
    super().__init__(probability)
  


  @classmethod
  def get_name(cls) -> str:
    return 'GA_Mutation_SingleSwitch'



  def mutate_sequence(self, sequence: str) -> str:
    """
    Muta la secuencia especificada por `sequence` y retorna la nueva secuencia.
    """
    i = random.randint(0, len(sequence) - 1)
    letter = sequence[i]
    while sequence[i] == letter:
      letter = random.choices(AMINOACIDS, k=1)[0]
    # en Python, las cadenas son inmutables
    temp = list(sequence)
    temp[i] = letter
    return ''.join(temp)
