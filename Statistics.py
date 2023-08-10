from dataclasses import dataclass, field, asdict
from .Population import Population
import statistics
from .Sequence import AMINOACIDS





@dataclass
class Statistics:

  iteration_id: int = field(default=None)
  min_fitness: float = field(default=None)
  fitness_mean: float = field(default=None)
  max_fitness: float = field(default=None)
  string_diversity: float = field(default=None)



  @classmethod
  def new_from_population(cls, population: Population):
    return cls(population.iterationId,
               min(population).fitness,
               statistics.fmean(population),
               max(population).fitness,
               cls.compute_string_diversity(population))
  


  def as_dict(self) -> dict:
    return asdict(self)



  @classmethod
  def compute_string_diversity(cls, population: Population) -> float:
    counts = []
    for i in range(len(population[0])):
      uniques = { individual[i] for individual in population }
      counts.append(len(uniques))
    return statistics.fmean(counts)
