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

  # TODO Revisar por que no coincide el optimo global con la secuencia objetivo

  _blosum_matrix = bl.BLOSUM(62)



  def __init__(self, 
               targetSequence: str,
               wingLength: int = 1,
               separation: int = 0) -> None:
    super().__init__({})

    if wingLength < 1: wingLength = 1
    if wingLength > 4: wingLength = 4
    if separation < 0: separation = 0
    if separation > 2: separation = 2

    self._target_sequence = targetSequence
    self._wing_length = wingLength
    self._separation = separation
    self._window_size = 2 * wingLength + 1
    self._sequence_length = len(targetSequence)
    self._residue_ordinals = self._compute_residue_ordinals()
    self._residue_window_indices = self._compute_residue_window_indices(targetSequence)
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
      'wingLength': self._wing_length,
      'separation': self._separation
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
    # para convertir una secuencia de aminoácidos a un vector de la misma
    # longitud que contiene números en el intervalo [-5.12, 5.12], vamos a 
    # apoyarnos de la matriz BLOSUM. El objetivo es lograr que, entre más se 
    # parezca una secuencia a la secuencia objetivo, mayor sea su aptitud. 
    # Entonces, por cada posición en la secuencia, vamos a ordenar los 
    # aminoácidos de manera ascendente según su score BLOSUM62, usando
    # como referencia el aminoácido correspondiente en la secuencia objetivo 
    residue_ordinals = []
    for res in self._target_sequence:
      amino_acid_scores = [
        (key, __class__._blosum_matrix[res][key])
        for key in __class__._blosum_matrix[res].keys()
        if key in AMINO_ACIDS_SET
      ]
      amino_acid_scores.sort(key=operator.itemgetter(1))
      residue_ordinals.append({
        item[0]: i
        for i, item in enumerate(amino_acid_scores)
      })
    return residue_ordinals
  


  def _compute_residue_window_indices(self, sequence: str) -> List[int]:
    # Al convertir un residue a un número real en el intervalo [-5.12, 5.12],
    # inevitablemente habrá una reducción de dimensionalidad, ya que el número
    # de posibles aminoácidos en cada residuo es finito. Este problema se 
    # amortigua al usar una subsecuencia o ventana (en lugar de un solo 
    # aminoácido) por cada residuo en la secuencia. Los índices de estas 
    # ventanas se precalculan con esta función
    distance = self._separation + 1
    residue_window_indices = []
    for pivot in range(len(sequence)):
      left_indices = [
        int((pivot - (distance * i)) % len(sequence))
        for i in reversed(range(1, self._wing_length + 1))
      ]
      right_indices = [
        int((pivot + (distance * i)) % len(sequence))
        for i in range(1, self._wing_length + 1)
      ]
      residue_window_indices.append(left_indices + [ pivot ] + right_indices)
    return residue_window_indices



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
    ind = self._residue_window_indices[pivot]
    return sum([
      self._residue_ordinals[ind[i]][sequence[ind[i]]] * 20**(len(ind) - 1 - i)
      for i in range(len(ind))
    ])

  

  def _to_real(self, decimal_residue: int) -> float:
    if decimal_residue % 2 == 0:
      return decimal_residue / 2 * self._STEP_SIZE
    return (decimal_residue + 1) / 2 * -self._STEP_SIZE
