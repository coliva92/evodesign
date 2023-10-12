from .Selection import Selection
from typing import List
from evodesign.Individual import Individual
from evodesign.Population import Population
import random
import evodesign.Choice as Choice





class Overselection(Selection):

  @classmethod
  def name(cls) -> str:
    return 'GA_Selection_Overselection'
  


  def __init__(self, 
               selectionSize: int,
               topSize: int,
               topTopProbability: float = 0.8,
               topBottomProbability: float = 0.0,
               bottomBottomProbability: float = 0.2) -> None:
    super().__init__(selectionSize)
    self._top_size = topSize
    self._top_top_probability = topTopProbability
    self._top_bot_probability = topBottomProbability
    self._bot_bot_probability = bottomBottomProbability
    self._weights = [
      topTopProbability,
      topBottomProbability,
      bottomBottomProbability
    ]
  


  def params_json(self) -> dict:
    params = super().params_json()
    params['topSize'] = self._top_size
    params['topTopProbability'] = self._top_top_probability
    params['topBottomProbability'] = self._top_bot_probability
    params['bottomBottomProbability'] = self._bot_bot_probability
    return params
  


  def select_parents(self, population: Population) -> List[Individual]:
    selected_parents = []
    top = population[-self._top_size:]
    bottom = population[:-self._top_size]
    while len(selected_parents) < self._selection_size:
      option = random.choices([ 0, 1, 2 ], self._weights)[0]
      if option == 0:
        mother, father = random.choices(top, k=2)
      if option == 1:
        mother = random.choice(top)
        father = random.choice(bottom)
      if option == 2:
        mother, father = random.choices(bottom, k=2)
      while mother.sequence == father.sequence:
        father = random.choice(top) if option == 0 else random.choice(bottom)
      selected_parents.append(mother)
      selected_parents.append(father)
    return selected_parents
