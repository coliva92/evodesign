from .FitnessFunction import FitnessFunction
from typing import Dict, List
from ..Sequence import _TEST_SEQUENCE, AMINOACIDS
import math
import random





class NegativeRastrigin(FitnessFunction):
  """[Función de Rastrigin](https://www.sfu.ca/~ssurjano/rastr.html) para 
  realizar pruebas sin utilizar directamente el predictor de la estructura de 
  una proteína al calcular la aptitud. Esta función debe utilizarse en conjunto 
  con `Predictor_Null`
  """
  
  _domain_values = [ 
    0.0, 0.5688, 1.1376, 1.7064, 2.2752, 2.844, 3.4128, 3.9816, 4.5504, 5.12, 
    -5.12, -4.608, -4.096, -3.584, -3.072, -2.56, -2.048, -1.536, -1.024, 
    -0.512 
  ]
  _target_sequence = None
  _aminoacid_ordinals = None



  def __init__(self) -> None:
    super().__init__({})



  @classmethod
  def get_name(cls) -> str:
    return 'Fitness_NegativeRastrigin'
  


  @classmethod
  def upper_bound(cls) -> float:
    return 0.0
  


  def compute_fitness(self, 
                      sequence: str, 
                      metrics: Dict[str, float]
                      ) -> float:
    x = self._to_rastrigin_domain(sequence)
    Sigma = 0.0
    for x_i in x:
      Sigma += (x_i ** 2 - 10.0 * math.cos(2.0 * math.pi * x_i))
    return -10.0 * len(sequence) - Sigma # queremos el valor negativo



  def _to_rastrigin_domain(self, sequence: str) -> List[float]:
    # La función de Rastrigin toma como entrada un vector x de
    # d componentes, donde cada componente x_i es un número real en el 
    # intervalo cerrado [-5.12, 5.12]. El óptimo global se encuentra en 
    # x = (0, ..., 0), y su valor es f(x) = 0.
    # El objetivo de `_to_rastrigin_domain` es transformar una secuencia de 
    # aminoácidos a un vector de números reales para la función de Rastrigin, 
    # de tal manera que el vector para el óptimo global corresponda con la 
    # secuencia de la estructura objetivo.
    # El problema es que la secuencia verdadera es desconocida para el 
    # algoritmo evolutivo. Así, vamos a ignorar la secuencia verdadera, y vamos 
    # a forzar el algoritmo a buscar una secuencia particular predeterminada.
    # Esto no debería afectar las pruebas con el algoritmo puesto que, de 
    # cualquier forma, ya se está ignorando la estructura objetivo al utilizar 
    # la función de prueba de Rastrigin. 
    if self._target_sequence == None:
      d = len(sequence)
      m = len(_TEST_SEQUENCE)
      temp = _TEST_SEQUENCE
      while d > m: temp += _TEST_SEQUENCE
      self._target_sequence = temp[:d]
      self._aminoacid_ordinals = []
      for _ in range(d):
        perm = random.sample(list(range(20)), 20)
        temp = {}
        for i, aa in enumerate(AMINOACIDS): temp[aa] = perm[i]
        self._aminoacid_ordinals.append(temp)
    x, n = [], len(sequence)
    for i in range(n):
      a = self._aminoacid_ordinals[i][self._target_sequence[i]]
      b = self._aminoacid_ordinals[i][sequence[i]]
      j = b - a
      if j > 9: j -= 20
      x.append(self._domain_values[j % 20])
    return x
