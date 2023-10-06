from .FitnessFunction import FitnessFunction
from ..Sequence import AMINO_ACIDS_SET
import blosum as bl
from typing import Dict, List
import math
import operator





class NegativeRastrigin(FitnessFunction):
  """[Función de Rastrigin](https://www.sfu.ca/~ssurjano/rastr.html) para 
  realizar pruebas sin utilizar directamente el predictor de la estructura de 
  una proteína al calcular la aptitud. Esta función debe utilizarse en conjunto 
  con `Predictor_Null`
  """

  _blosum_matrix = bl.BLOSUM(62)
  _domain_values = [ 
    -5.12, 5.12, -4.608, 4.5504, -4.096, 3.9816, -3.584, 3.4128, -3.072, 2.844,
    -2.56, 2.2752, -2.048, 1.7064, -1.536, 1.1376, -1.024, 0.5688, -0.512,
    0.0
  ]



  def __init__(self, targetSequence: str) -> None:
    super().__init__({})
    self._target_sequence = targetSequence
    self._sequence_length = len(targetSequence)
    self._residue_scores = self._compute_residue_scores()



  @classmethod
  def get_name(cls) -> str:
    return 'Fitness_NegativeRastrigin'



  @classmethod
  def upper_bound(cls) -> float:
    return 0.0



  def compute_fitness(self, 
                      sequence: str, 
                      _: Dict[str, float]
                      ) -> float:
    x = self._to_rastrigin_domain(sequence)
    sigma = sum([ 
      x_i ** 2 - 10.0 * math.cos(2.0 * math.pi * x_i) 
      for x_i in x 
    ])
    return -10.0 * self._sequence_length - sigma # queremos el valor negativo
  


  def _compute_residue_scores(self) -> List[Dict[str, float]]:
    residue_scores = []
    for res in self._target_sequence:
      temp = [
        (key, __class__._blosum_matrix[res][key])
        for key in __class__._blosum_matrix[res].keys()
        if key in AMINO_ACIDS_SET
      ]
      temp = sorted(temp, key=operator.itemgetter(1))
      residue_scores.append({
        item[0]: __class__._domain_values[i]
        for i, item in enumerate(temp)
      })
    return residue_scores



  def _to_rastrigin_domain(self, sequence: str) -> List[float]:
    # La función de Rastrigin toma como entrada un vector x de
    # d componentes, donde cada componente x_i es un número real en el 
    # intervalo cerrado [-5.12, 5.12]. El óptimo global se encuentra en 
    # x = (0, ..., 0), y su valor es f(x) = 0.
    # El objetivo de `_to_rastrigin_domain` es transformar una secuencia de 
    # aminoácidos a un vector de números reales para la función de Rastrigin, 
    # de tal manera que el vector para el óptimo global corresponda con la 
    # secuencia de la estructura objetivo.
    return [ self._residue_scores[i][res] for i, res in enumerate(sequence) ]
