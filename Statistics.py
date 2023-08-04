from dataclasses import dataclass, field, asdict
from .Population import Population
import statistics





@dataclass
class Statistics:

  iteration_id: int = field(default=None)
  min_fitness: float = field(default=None)
  fitness_mean: float = field(default=None)
  max_fitness: float = field(default=None)



  @classmethod
  def new_from_population(cls, population: Population):
    return cls(population.iterationId,
               min(population).fitness,
               statistics.fmean(population),
               max(population).fitness)
  


  def as_dict(self) -> dict:
    return asdict(self)
