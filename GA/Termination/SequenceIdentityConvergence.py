from .Terminator import Terminator
from evodesign.Statistics import Statistics
from evodesign.Population import Population
import math





class SequenceIdentityConvergence(Terminator):
  
  def __call__(self, 
               population: Population,
               stats: Statistics
               ) -> bool:
    if not stats:
      stats = Statistics.new_from_population(population)
    n = len(population[0])
    return stats.sequence_identity >= math.floor(0.9 * n)
