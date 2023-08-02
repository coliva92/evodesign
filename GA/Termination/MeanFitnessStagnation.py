from .Terminator import Terminator
from typing import List, Optional
from evodesign.Statistics import Statistics
from evodesign.Individual import Individual
import os
import csv





class MeanFitnessStagnation(Terminator):
  
  def __init__(self, statsFilename: Optional[str] = None) -> None:
    super().__init__()
    self._avg_fitnesses = []
    if statsFilename and os.path.isfile(statsFilename):
      with open(statsFilename, 'rt', encoding='utf-8') as csv_file:
        stats = csv.reader(csv_file, delimiter=',', quotechar='"')
        next(stats)
        self._avg_fitnesses = [ row[2] for row in stats ]
  


  def __call__(self, 
               iterationId: int, 
               population: List[Individual],
               stats: Statistics
               ) -> bool:
    self._avg_fitnesses.append(stats.fitness_mean)
    if iterationId >= 10:
      older_avg = self._avg_fitnesses[iterationId - 10]
      newer_avg = self._avg_fitnesses[iterationId]
      return newer_avg - older_avg < 0.0001
    return False
