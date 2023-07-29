from .Selection import Selection
from typing import List
from evodesign.Population import Individual
import random





class Overselection(Selection):
  """
  Operación de selección donde los primeros `topSize` individuos de la 
  población se seleccionan como padres con una probabilidad de 
  `topProbability`, y donde los demás individuos se seleccionan con una 
  probabilidad de `1.0 - topProbability`.
  """

  _options = [ True, False ]



  @classmethod
  def get_name(cls) -> str:
    return 'GA_Selection_Overselection'
  


  def __init__(self, 
               selectionSize: int,
               topSize: int,
               topProbability: float = 0.8) -> None:
    """
    Constructor.
    - `selectionSize`: el número de individuos a elegir como padres.
    - `topSize`: el tamaño del grupo de los mejores _N_ individuos; la 
      probabilidad de selección está sesgada hacia los individuos de este grupo.
    - `topProbability`: la probabilidad de elegir como padre un individuo del 
      grupo de los primeros _N_ individuos. 
    """
    super().__init__(selectionSize)
    self._top_size = topSize
    self._weights = [ topProbability, 1.0 - topProbability ]
    self._top_probability = topProbability
  


  def get_params_memento(self) -> dict:
    params = super().get_params_memento()
    params['topSize'] = self._top_size
    params['topProbability'] = self._top_probability
    return params
  


  def select_parents(self, population: List[Individual]) -> List[Individual]:
    """
    Selecciona algunos individuos del arreglo especificado por `population` 
    para posteriormente recombinarlos y generar la población de la siguiente 
    generación.
    """
    if random.choices(Overselection._options, self._weights, k=1)[0]:
      return random.sample(population[-self._top_size:], self._selection_size)
    return random.sample(population[0:-self._top_size], self._selection_size)
