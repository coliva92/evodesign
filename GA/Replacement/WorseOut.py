from .Replacement import Replacement
from typing import List
from evodesign.Population import Population





class WorseOut(Replacement):
  
  @classmethod
  def name(cls) -> str:
    return 'GA_Replacement_WorseOut'



  def __call__(self, 
               population: Population,
               children: Population
               ) -> Population:
    # insertamos los hijos a la población, preservando el orden del segundo
    next_population = population.individuals.copy()
    for child in children:
      for i, individual in enumerate(next_population):
        if child <= individual:
          next_population.insert(i, child)
          break
    return Population(next_population[len(children):],
                      population.iteration_id + 1)
