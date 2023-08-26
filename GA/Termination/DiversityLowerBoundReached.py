from .Terminator import Terminator
from evodesign.Statistics import Statistics
from evodesign.Population import Population
import math





class DiversityLowerBoundReached(Terminator):
  
  def __call__(self, 
               population: Population,
               stats: Statistics
               ) -> bool:
    return stats.diversity < math.ceil(0.1 * len(population[0]))
