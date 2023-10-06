from .Terminator import Terminator
from evodesign.Statistics import Statistics
from evodesign.Population import Population
import math





class DiversityLowerBoundReached(Terminator):
  
  def __call__(self, 
               population: Population,
               stats: Statistics
               ) -> bool:
    if not stats:
      stats = Statistics.new_from_population(population)
    return stats.sequence_diversity < math.floor(0.1 * len(population[0]))
