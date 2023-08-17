from dataclasses import dataclass, field, asdict
from typing import Tuple
from .Population import Population
import statistics
import math





@dataclass
class Statistics:

  iteration_id: int = field(default=None)
  min_fitness: float = field(default=None)
  fitness_mean: float = field(default=None)
  max_fitness: float = field(default=None)
  diversity: float = field(default=None)



  @classmethod
  def new_from_population(cls, population: Population):
    smallest, largest, average = cls.min_max_mean(population)
    diversity = cls.average_hamming_distance(population)
    return cls(population.iteration_id,
               smallest,
               average,
               largest,
               diversity)
  


  @classmethod
  def min_max_mean(cls, population: Population) -> Tuple[float, float, float]:
    # una prueba empírica con una población de 1M de secuencias mostró que usar
    # las funciones min, max y statistics.fmean es aprox. siete veces más lento
    # que calcular estos valores en un solo bucle
    summation, smallest, largest = 0.0, math.inf, -math.inf
    for individual in population:
      if individual.fitness < smallest:
        smallest = individual.fitness
      if individual.fitness > largest:
        largest = individual.fitness
      summation += individual.fitness
    average = summation / len(population)
    return smallest, largest, average
  


  @classmethod
  def average_hamming_distance(cls, population: Population) -> float:
    # una prueba con una población de 1K secuencias de longitud 1K mostró que 
    # no hay diferencia relevante en tiempo de ejecución entre usar list 
    # comprehensions y usar bucles for; usamos bucles for para hacer el código 
    # más legible
    averages = []
    for i, current in enumerate(population[:-1]):
      counts = []
      for other in population[i + 1:]:
        counts.append(sum(map(lambda pair: pair[0] != pair[1], 
                              zip(current, other))))
      averages.append(statistics.fmean(counts) / len(population[:-1]))
    return sum(averages)
  


  def as_dict(self) -> dict:
    return asdict(self)
