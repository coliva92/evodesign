from .Metric import Metric
from typing import Dict, List, Optional
import evodesign.Sequence as Sequence
import evodesign.Chain as Chain
from Bio.Align import substitution_matrices
import math
import operator
import pandas as pd
from Bio.PDB.Atom import Atom
from ..Context import Context





class Rastrigin(Metric):

  """[Función de Rastrigin](https://www.sfu.ca/~ssurjano/rastr.html) para 
  realizar pruebas sin utilizar directamente el predictor de la estructura de 
  una proteína al calcular la aptitud. Esta función debe utilizarse en conjunto 
  con `Predictor_Null`
  """
  


  _VALID_WINDOW_WIDTH_VALUES = { 1, 3, 5, 7, 9 }
  def _params(self) -> dict:
    params = super()._params()
    params['window_width'] = self._window_width
    return params
  


  def __init__(self,
               column: Optional[str] = None
               ) -> None:
    super().__init__(column)
    self._window_width = 3
    self._wing_size = (self._window_width - 1) // 2
    self._residue_ordinals = None
    self._STEP_SIZE = 2 * 5.12 / 20 ** self._window_width



  def _compute_values(self, 
                      backbone: List[Atom],
                      data: pd.Series,
                      context: Context
                      ) -> pd.Series:
    if not self._residue_ordinals:
      self._target_sequence = context.ref_sequence
      blosum_matrix = substitution_matrices.load("BLOSUM62")
      self._residue_ordinals = self._compute_residue_ordinals(blosum_matrix)
    x = self._to_rastrigin_domain(data['sequence'])
    sigma = sum([ 
      x_i ** 2 - 10.0 * math.cos(2.0 * math.pi * x_i) 
      for x_i in x 
    ])
    data[self.column_name()] = -10.0 * len(data['sequence']) - sigma # queremos el valor negativo
    return data
  


  def _compute_residue_ordinals(self, 
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
    for residue in self._target_sequence:
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
    decimals = [ self._to_decimal(sequence, i) for i in range(len(sequence)) ]
    reals = [ self._to_real(d) for d in decimals ]
    return reals
  


  def _to_decimal(self, 
                  sequence: str, 
                  pivot: int
                  ) -> int:
    # Dado que el conjunto de aminoácidos que se permiten en cada residuo es
    # discreto y finito, si deseamos asociar cada aminoácido en dicho conjunto 
    # con un número real en el intervalo [-5.12, 5.12], inevitablemente algunos 
    # de los números en este intervalo no podrán ser representados (por la 
    # discretización). Este problema se amortigua al usar una subsecuencia o 
    # ventana (en lugar de un solo aminoácido) por cada posición o residuo en 
    # la secuencia.
    n = len(sequence)
    a = pivot - self._wing_size
    b = pivot + self._wing_size + 1
    indices = [ (a - 5) % n ] + \
              [ (a - 2) % n ] + \
              [ i % n for i in range(a, b) ] + \
              [ (b + 2) % n ] + \
              [ (b + 5) % n ]
    ordinals = [
      self._residue_ordinals[j][sequence[j]] * 20**(self._window_width - 1 - i)
      for i, j in enumerate(indices)
    ]
    return sum(ordinals)

  

  def _to_real(self, decimal_residue: int) -> float:
    if decimal_residue % 2 == 0:
      return decimal_residue / 2 * self._STEP_SIZE
    return (decimal_residue + 1) / 2 * -self._STEP_SIZE
