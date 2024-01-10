from .Terminator import Terminator
from typing import Optional
from evodesign.Statistics import Statistics
from evodesign.Population import Population
import os
import csv





class FitnessMeanConvergence(Terminator):
  
  def __init__(self, statsFilename: Optional[str] = None) -> None:
    super().__init__()
    self._avg_fitnesses = []
    if statsFilename and os.path.isfile(statsFilename):
      with open(statsFilename, 'rt', encoding='utf-8') as csv_file:
        stats = csv.DictReader(csv_file, dialect='unix')
        self._avg_fitnesses = [ float(row['fitness_mean']) for row in stats ]
  


  def __call__(self, 
               population: Population,
               stats: Statistics
               ) -> bool:
    self._avg_fitnesses.append(stats.fitness_mean)
    if population.iteration_id >= 10:
      older_avg = self._avg_fitnesses[population.iteration_id - 10]
      newer_avg = self._avg_fitnesses[population.iteration_id]
      return newer_avg - older_avg < 0.0001
    return False
