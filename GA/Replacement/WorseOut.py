from .Replacement import Replacement
from typing import List
from evodesign.Population import Population
from evodesign.Individual import Individual





class WorseOut(Replacement):
  
  @classmethod
  def get_name(cls) -> str:
    return 'GA_Replacement_WorseOut'



  def __call__(self, 
               population: Population,
               children: List[Individual]
               ) -> Population:
    next_population = population.individuals.copy()
    for child in children:
      for i, individual in enumerate(next_population):
        if child <= individual:
          next_population.insert(i, child)
          break
    return Population(population.iterationId + 1, 
                      next_population[len(children):])
