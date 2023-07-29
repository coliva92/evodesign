from .Mutation import Mutation
import random
from evodesign.Sequence import AMINOACIDS





class MultipleSwitches(Mutation):
  """
  Operación de mutación donde se elige aleatoriamente un residuo en la 
  secuencia y se cambia por otro aminoácido aleatorio.
  """

  @classmethod
  def get_name(cls) -> str:
    return 'GA_Mutation_MultipleSwitches'
  
  

  def __init__(self, probability: float = 1.0, numSwitches: int = 1) -> None:
    """
    Constructor.
    - `probability`: la probabilidad con la que se aplicará la mutación sobre 
      un individuo.
    """
    super().__init__(probability)
    self._num_switches = numSwitches
  


  def get_params_memento(self) -> dict:
    params = super().get_params_memento()
    params['numSwitches'] = self._num_switches
    return params



  def mutate_sequence(self, sequence: str) -> str:
    """
    Muta la secuencia especificada por `sequence` y retorna la nueva secuencia.
    """
    positions = random.sample(list(range(len(sequence))), self._num_switches)
    # en Python, las cadenas son inmutables; hay que convertirlas a una lista
    seq_list = list(sequence)
    for i in positions:
      letter = seq_list[i]
      while seq_list[i] == letter: 
        letter = random.choice(AMINOACIDS)
      seq_list[i] = letter
    return ''.join(seq_list)
