from . import Terminator
from ... import Statistics
from ... import Population
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
