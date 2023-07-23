from .Mutation import Mutation
import random
from evodesign.Sequence import AMINOACID_CATEGORIES, AMINOACIDS_BY_CATEGORY





class ChemicalSwitch(Mutation):
  """
  Operación de mutación donde se elige aleatoriamente un residuo en la 
  secuencia y se cambia por otro aminoácido aleatorio.
  """

  _options = [ True, False ]
  _category_complements = {
    'negative': [ 'positive', 'unchargedPolar', 'nonpolar' ],
    'positive': [ 'negative', 'unchargedPolar', 'nonpolar' ],
    'unchargedPolar': [ 'negative', 'positive', 'nonpolar' ],
    'nonpolar': [ 'negative', 'positive', 'unchargedPolar' ]
  }



  def __init__(self, 
               probability: float = 1.0, 
               preservationBias: float = 0.5) -> None:
    """
    Constructor.
    - `probability`: la probabilidad con la que se aplicará la mutación sobre 
      un individuo.
    """
    super().__init__(probability)
    self._weights = [ preservationBias, 1.0 - preservationBias ]
    self._preservation_bias = preservationBias



  @classmethod
  def get_name(cls) -> str:
    return 'GA_Mutation_ChemicalSwitch'
  


  def get_params_memento(self) -> dict:
    params = super().get_params_memento()
    params['preservationBias'] = self._preservation_bias
    return params



  def mutate_sequence(self, sequence: str) -> str:
    """
    Muta la secuencia especificada por `sequence` y retorna la nueva secuencia.
    """
    i = random.randint(0, len(sequence) - 1)
    letter = sequence[i]
    preserve_category = random.choices(self._options, self._weights, k=1)[0]
    if preserve_category:
      options = AMINOACIDS_BY_CATEGORY[AMINOACID_CATEGORIES[letter]]
      while sequence[i] == letter:
        letter = random.choice(options)
    else:
      options = self._category_complements[AMINOACID_CATEGORIES[letter]]
      category = random.choice(options)
      letter = random.choice(AMINOACIDS_BY_CATEGORY[category])
    # en Python, las cadenas son inmutables
    temp = list(sequence)
    temp[i] = letter
    return ''.join(temp)
