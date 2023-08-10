from dataclasses import dataclass, field, asdict
from .Population import Population
import statistics
from Sequence import AMINOACIDS





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
               cls._compute_string_diversity(population))
  


  def as_dict(self) -> dict:
    return asdict(self)



  @classmethod
  def _compute_string_diversity(cls, population: Population) -> float:
    n = len(AMINOACIDS) if len(population) >= len(AMINOACIDS) \
                        else len(population)
    ratios = []
    for i in range(len(population[0])):
      uniques = { individual[i] for individual in population }
      ratios.append(len(uniques) / n)
    return statistics.fmean(ratios)
