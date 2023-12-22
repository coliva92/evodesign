from .Terminator import Terminator
from evodesign.Statistics import Statistics
from evodesign.Population import Population





class SequenceIdentityConvergence(Terminator):
  
  def __call__(self, 
               population: Population,
               stats: Statistics
               ) -> bool:
    if not stats:
      stats = Statistics.new_from_population(population)
    return stats.sequence_identity > 0.9
