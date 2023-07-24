from .Mutation import Mutation
import random
from evodesign.Sequence import AMINOACIDS





class RandomResetting(Mutation):
  """
  Operación de mutación donde se elige aleatoriamente un residuo en la 
  secuencia y se cambia por otro aminoácido aleatorio.
  """

  def __init__(self, 
               probability: float = 1.0, 
               residueProbability: float = 0.1) -> None:
    """
    Constructor.
    - `probability`: la probabilidad con la que se aplicará la mutación sobre 
      un individuo.
    """
    super().__init__(probability)
    self._residue_prob = residueProbability
    self._weights = [ residueProbability, 1.0 - residueProbability ]
  


  @classmethod
  def get_name(cls) -> str:
    return 'GA_Mutation_RandomResetting'
  


  def get_params_memento(self) -> dict:
    params = super().get_params_memento()
    params['residueProbability'] = self._residue_prob
    return params



  def mutate_sequence(self, sequence: str) -> str:
    """
    Muta la secuencia especificada por `sequence` y retorna la nueva secuencia.
    """
    temp = list(sequence)
    for i in range(len(temp)):
      if random.choices(Mutation._options, self._weights)[0]:
        letter = temp[i]
        while temp[i] == letter:
          letter = random.choice(AMINOACIDS)
        temp[i] = letter
    return ''.join(temp)
