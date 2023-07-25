from typing import Dict
from typing import Optional
import evodesign.Sequence as Sequence





class Individual:
  """
  La representación de un individuo (i.e., una secuencia de aminoácidos) en la 
  población del algoritmo evolutivo.
  """

  _next_id = 0



  def __init__(self, 
               sequence: str, 
               id: Optional[str] = None,
               fitness: Optional[float] = None,
               metrics: Optional[Dict[str, float]] = None) -> None:
    """
    Constructor.
    - `sequence`: la secuencia de aminoácidos que representa este individuo.
    - `id`: el nombre con el que se identificará este individuo en el algoritmo 
      evolutivo.
    - `fitness`: el valor de aptitud para este individuo.
    - `metrics`: un diccionario que contiene valores de diferentes métricas de 
      calidad para este individuo; estos valores se utilizan para calcular la 
      aptitud del individuo.
    """
    if metrics == None: metrics = {}
    self.id = id
    if id != None:
      integer_id = int(float(id))
      if integer_id >= Individual._next_id:
        Individual._next_id = integer_id + 1
    else:
      self.id = Individual._create_id()
    self.sequence = sequence
    self.fitness = fitness
    self.metrics = metrics
  


  @classmethod
  def random(cls, sequenceLength: int):
    sequence = Sequence.create_random_sequence(sequenceLength)
    return cls(sequence)
  


  def get_memento(self) -> dict:
    memento = { 
      'id': self.id,
      'fitness': self.fitness,
      'metrics': {},
      'sequence': self.sequence
    }
    for key, value in self.metrics.items():
      memento['metrics'][key] = value
    return memento
    

  
  def __eq__(self, other) -> bool:
    if not isinstance(other, Individual):
      return False
    return self.fitness == other.fitness
  


  def __lt__(self, other):
    if not isinstance(other, Individual):
      return True
    return self.fitness < other.fitness
  


  def __le__(self, other):
    if not isinstance(other, Individual):
      return True
    return (self.fitness < other.fitness) or (self.fitness ==  other.fitness)



  @classmethod
  def _create_id(cls) -> str:
    id = f'{cls._next_id}'
    cls._next_id += 1
    return id
