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



  def __init__(self, 
               targetSequence: str,
               windowSize: int = 3,
               distance: int = 1) -> None:
    super().__init__({})
    if windowSize % 2 == 0:
      raise RuntimeError('window size must be an odd integer')
    self._target_sequence = targetSequence
    self._window_size = windowSize
    self._distance = distance
    self._sequence_length = len(targetSequence)
    self._residue_ordinals = self._compute_residue_ordinals()
    self._STEP_SIZE = 2 * 5.12 / 20 ** self._window_size



  @classmethod
  def name(cls) -> str:
    return 'Fitness_NegativeRastrigin'



  @classmethod
  def upper_bound(cls) -> float:
    return 0.0
  


  def params_json(self) -> dict:
    return {
      'targetSequence': self._target_sequence,
      'windowSize': self._window_size,
      'distance': self._distance
    }



  def compute_fitness(self,
                      _: Dict[str, float],
                      sequence: str
                      ) -> float:
    x = self._to_rastrigin_domain(sequence)
    sigma = sum([ 
      x_i ** 2 - 10.0 * math.cos(2.0 * math.pi * x_i) 
      for x_i in x 
    ])
    return -10.0 * self._sequence_length - sigma # queremos el valor negativo
  


  def _compute_residue_ordinals(self) -> List[Dict[str, float]]:
    residue_ordinals = []
    for res in self._target_sequence:
      amino_acid_scores = [
        (key, __class__._blosum_matrix[res][key])
        for key in __class__._blosum_matrix[res].keys()
        if key in AMINO_ACIDS_SET
      ]
      amino_acid_scores = sorted(amino_acid_scores, key=operator.itemgetter(1))
      residue_ordinals.append({
        item[0]: i
        for i, item in enumerate(amino_acid_scores)
      })
    return residue_ordinals



  def _to_rastrigin_domain(self, sequence: str) -> List[float]:
    # La función de Rastrigin toma como entrada un vector x de
    # d componentes, donde cada componente x_i es un número real en el 
    # intervalo cerrado [-5.12, 5.12]. El óptimo global se encuentra en 
    # x = (0, ..., 0), y su valor es f(x) = 0.
    # El objetivo de `_to_rastrigin_domain` es transformar una secuencia de 
    # aminoácidos a un vector de números reales para la función de Rastrigin, 
    # de tal manera que el vector para el óptimo global corresponda con la 
    # secuencia de la estructura objetivo.
    return [ 
      self._to_real(self._to_decimal(sequence, i)) 
      for i in range(len(sequence)) 
    ]
  


  def _to_decimal(self, sequence: str, pivot: int) -> int:
    half_window = (self._window_size - 1) // 2
    indexes = []
    for i in reversed(range(1, half_window + 1)):
      idx = int((pivot - (self._distance * i)) % len(sequence))
      indexes.append(idx)
    indexes.append(pivot)
    for i in range(1, half_window + 1):
      idx = int((pivot + (self._distance * i)) % len(sequence))
      indexes.append(idx)
    return sum([
      self._residue_ordinals[indexes[i]][sequence[indexes[i]]] * 20 ** (len(indexes) - 1 - i)
      for i in range(len(indexes))
    ])
    # TODO borrar este codigo:
    r = (pivot - 1) % len(sequence)
    s = (pivot + 1) % len(sequence)
    return self._residue_ordinals[r][sequence[r]] * 20 ** 2 + \
      self._residue_ordinals[pivot][sequence[pivot]] * 20 ** 1 + \
      self._residue_ordinals[s][sequence[s]]
  


  def _to_real(self, decimal_residue: int) -> float:
    if decimal_residue % 2 == 0:
      return decimal_residue / 2 * self._STEP_SIZE
    return (decimal_residue + 1) / 2 * -self._STEP_SIZE
