from typing import Dict, Optional
import evodesign.Sequence as Sequence
import time





class Individual:
  """
  La representación de un individuo (i.e., una secuencia de aminoácidos) en la 
  población del algoritmo evolutivo.
  """

  _next_id = 0



  @classmethod
  def random(cls, sequenceLength: int):
    return cls(Sequence.create_random_sequence(sequenceLength))



  def __init__(self, 
               sequence: str,
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
    if metrics is None: metrics = {}
    self.sequence = sequence
    self.fitness = fitness
    self.metrics = metrics
  


  def get_memento(self) -> dict:
    memento = {
      'fitness': self.fitness,
      'metrics': {},
      'sequence': self.sequence
    }
    for key, value in self.metrics.items():
      memento['metrics'][key] = value
    return memento
    

  
  def __eq__(self, other) -> bool:
    if other.__class__ != self.__class__:
      return NotImplemented
    return self.fitness == other.fitness
  


  def __lt__(self, other):
    if other.__class__ != self.__class__:
      return NotImplemented
    return self.fitness < other.fitness
  


  def __le__(self, other):
    if other.__class__ != self.__class__:
      return NotImplemented
    return (self.fitness < other.fitness) or (self.fitness ==  other.fitness)
