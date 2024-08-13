from .Metric import Metric
from typing import Dict, List, Optional
import evodesign.Sequence as Sequence
from Bio.Align import substitution_matrices
import math
import operator
import pandas as pd
from Bio.PDB.Atom import Atom
from ..Context import Context
import numpy as np





class Rastrigin(Metric):

  """[Función de Rastrigin](https://www.sfu.ca/~ssurjano/rastr.html) para 
  realizar pruebas sin utilizar directamente el predictor de la estructura de 
  una proteína al calcular la aptitud. Esta función debe utilizarse en conjunto 
  con `Predictor_Null`
  """
  


  def __init__(self,
               column: Optional[str] = None
               ) -> None:
    super().__init__(column)
    self._STEP_SIZE = 2 * 5.12 / 20 ** 6
    self._residue_ordinals = None
    self._subsequence_indices = None
    self._subsequence_ordinals = None
    self._rng = np.random.default_rng()



  def _compute_values(self, 
                      backbone: List[Atom],
                      data: pd.Series,
                      context: Context
                      ) -> pd.Series:
    if self._subsequence_indices is None:
      n = len(data["sequence"])
      self._subsequence_indices = [
        self._get_subsequence_indices(i, n)
        for i in range(len(data["sequence"]))
      ]
    if self._residue_ordinals is None:
      blosum_matrix = substitution_matrices.load("BLOSUM62")
      self._residue_ordinals = self._compute_residue_ordinals(context.ref_sequence, 
                                                              blosum_matrix)
    if self._subsequence_ordinals is None:
      self._subsequence_ordinals = [
        [
          self._residue_ordinals[indices[0]],
          self._residue_ordinals[indices[1]],
          self._reverse_ordinals(self._residue_ordinals[indices[2]]),
          self._reverse_ordinals(self._residue_ordinals[indices[3]]),
          self._residue_ordinals[indices[4]],
          self._reverse_ordinals(self._residue_ordinals[indices[5]])
        ]
        for indices in self._subsequence_indices
      ]
    x = self._to_rastrigin_domain(data['sequence'])
    data[self.column_name()] = self._compute_rastrigin_value(len(data['sequence']), x)
    return data
  


  def _compute_rastrigin_value(self, d: int, x: List[float]) -> float:
    sigma = sum([ 
      x_i ** 2 - 10.0 * math.cos(2.0 * math.pi * x_i) 
      for x_i in x 
    ])
    return -10.0 * d - sigma # queremos el valor negativo



  def _compute_residue_ordinals(self, 
                                target_sequence: str,
                                blosum_matrix: substitution_matrices.Array
                                ) -> List[Dict[str, float]]:
    # Para convertir una secuencia de aminoácidos a un vector de la misma
    # longitud que contiene números en el intervalo [-5.12, 5.12], vamos a 
    # apoyarnos de la matriz BLOSUM. El objetivo es lograr que, entre más se 
    # parezca una secuencia a la secuencia objetivo, mayor sea su aptitud. 
    # Entonces, por cada posición en la secuencia, vamos a ordenar los 
    # aminoácidos de manera ascendente según su score BLOSUM62, usando
    # como referencia el aminoácido correspondiente en la secuencia objetivo. 
    residue_ordinals = []
    for residue in target_sequence:
      amino_acid_scores = [
        (amino_acid, blosum_matrix[residue][amino_acid])
        for amino_acid in blosum_matrix[residue].keys()
        if amino_acid in Sequence.AMINO_ACIDS
      ]
      amino_acid_scores.sort(key=operator.itemgetter(1), reverse=True)
      residue_ordinals.append({
        score[0]: i # la clave es el aminoácido
        for i, score in enumerate(amino_acid_scores)
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
    subsequences = [ self._get_subsequence(sequence, i) for i in range(len(sequence)) ]
    decimals = [ self._to_decimal(i, subseq) for i, subseq in enumerate(subsequences) ]
    reals = [ self._to_real(d) for d in decimals ]
    return reals
  


  def _get_subsequence(self, 
                       sequence: str, 
                       idx: int
                       ) -> str:
    subsequence = [ sequence[i] for i in self._subsequence_indices[idx] ]
    return ''.join(subsequence)



  def _to_decimal(self, idx: int, sequence: str) -> int:
    # Dado que el conjunto de aminoácidos que se permiten en cada residuo es
    # discreto y finito, si deseamos asociar cada aminoácido en dicho conjunto 
    # con un número real en el intervalo [-5.12, 5.12], inevitablemente algunos 
    # de los números en este intervalo no podrán ser representados (por la 
    # discretización). Este problema se amortigua al usar una subsecuencia o 
    # ventana (en lugar de un solo aminoácido) por cada posición o residuo en 
    # la secuencia.
    digits = [
      self._subsequence_ordinals[idx][j][sequence[j]] * 20 ** (5 - j)
      for j in range(len(sequence))
    ]
    return sum(digits)

  

  def _to_real(self, decimal_residue: int) -> float:
    if decimal_residue % 2 == 0:
      return decimal_residue / 2 * self._STEP_SIZE
    return (decimal_residue + 1) / 2 * -self._STEP_SIZE



  def _get_subsequence_indices(self, idx: int, N: int) -> List[int]:
    neighbor_indices = [ idx, (idx - 1) % N, (idx + 1) % N ]
    indices = list(range(N))
    long_distance_indices = self._rng.choice(indices, size=3, replace=False)
    for i in range(3):
      while long_distance_indices[i] in neighbor_indices:
        long_distance_indices[i] = self._rng.choice(indices)
    return [
      neighbor_indices[0], neighbor_indices[1], long_distance_indices[0],
      long_distance_indices[1], neighbor_indices[2], long_distance_indices[2]
    ]
  


  def _shuffle_ordinals(self, 
                        ordinals: Dict[str, float]
                        ) -> Dict[str, float]:
    indices = list(range(1, 20))
    sorted_amino_acids = { i: letter for letter, i in ordinals.items() }
    shuffled_positions = self._rng.choice(indices, size=19, replace=False)
    amino_acids = [ sorted_amino_acids[0] ] + \
                  [ sorted_amino_acids[i] for i in shuffled_positions ]
    return { letter: i for i, letter in enumerate(amino_acids) }
  


  def _reverse_ordinals(self, ordinals: Dict[str, float]) -> Dict[str, float]:
    indices = list(reversed(range(1, 20)))
    sorted_amino_acids = { i: letter for letter, i in ordinals.items() }
    reversed_amino_acids = [ sorted_amino_acids[0] ] + \
                           [ sorted_amino_acids[i] for i in indices ]
    return { letter: i for i, letter in enumerate(reversed_amino_acids) }
