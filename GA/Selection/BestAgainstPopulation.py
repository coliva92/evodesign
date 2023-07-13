from .Selection import Selection
from typing import List
from evodesign.Population import Individual
import random





class BestAgainstPopulation(Selection):
  
  def __init__(self, selectionSize: int = 0) -> None:
    """
    Constructor.
    - `selectionSize`: el número de individuos a elegir como padres.
    - `topSize`: el tamaño del grupo de los mejores _N_ individuos; la 
      probabilidad de selección está sesgada hacia los individuos de este grupo.
    - `topProbability`: la probabilidad de elegir como padre un individuo del 
      grupo de los primeros _N_ individuos. 
    """
    super().__init__(0) # selectionSize no se utiliza
  


  @classmethod
  def get_name(cls) -> str:
    return 'GA_Selection_BestAgainstPopulation'
  


  def select_parents(self, population: List[Individual]) -> List[Individual]:
    """
    Selecciona algunos individuos del arreglo especificado por `population` 
    para posteriormente recombinarlos y generar la población de la siguiente 
    generación.
    """
    n = len(population) - 1
    best = population[-1]
    everyone_else = population[:-1]
    parents = []
    for i in range(n):
      parents.append(best)
      parents.append(everyone_else[i])
    return parents
