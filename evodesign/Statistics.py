from dataclasses import dataclass, field, asdict
from typing import Tuple
from Population import Population
import statistics
import math





@dataclass
class Statistics:

  iteration_id: int = field(default=None)
  min_fitness: float = field(default=None)
  fitness_mean: float = field(default=None)
  max_fitness: float = field(default=None)
  sequence_identity: float = field(default=None)
  residue_identity: float = field(default=None)



  @classmethod
  def new_from_population(cls, population: Population):
    smallest, largest, average = cls.min_max_mean(population)
    seq_diversity = cls.average_sequence_identity(population)
    res_diversity = cls.average_per_residue_identity(population)
    return cls(population.iteration_id,
               smallest,
               average,
               largest,
               seq_diversity,
               res_diversity)
  


  @classmethod
  def min_max_mean(cls, population: Population) -> Tuple[float, float, float]:
    # una prueba empírica con una población de 1M de secuencias mostró que usar
    # las funciones min, max y statistics.mean es aprox. siete veces más lento
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
  def average_per_residue_identity(cls, population: Population) -> float:
    return statistics.mean([
      20 - len({ sequence[i] for sequence in population })
      for i in range(len(population[0]))
    ])
  


  @classmethod
  def average_sequence_identity(cls, population: Population) -> float:
    return statistics.mean([
      statistics.mean([
        sum(map(lambda pair: pair[0] == pair[1], zip(current, other)))
        for other in population[i + 1:]
      ])
      for i, current in enumerate(population[:-1])
    ])
  


  def as_dict(self) -> dict:
    return asdict(self)
