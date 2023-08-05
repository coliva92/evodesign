from .Mutation import Mutation
import random
from evodesign.Sequence import AMINOACIDS
import evodesign.Choice as Choice





class RandomResetting(Mutation):

  @classmethod
  def get_name(cls) -> str:
    return 'GA_Mutation_RandomResetting'
  
  

  def __init__(self, 
               probability: float = 1.0, 
               residueProbability: float = 0.1
               ) -> None:
    super().__init__(probability)
    self._residue_prob = residueProbability
    self._weights = [ residueProbability, 1.0 - residueProbability ]
  


  def as_json(self) -> dict:
    params = super().as_json()
    params['residueProbability'] = self._residue_prob
    return params



  def mutate_sequence(self, sequence: str) -> str:
    def switch(old_letter: str) -> str:
      if not Choice.flip_coin(self._weights):
        return old_letter
      new_letter = random.choice(AMINOACIDS)
      while new_letter == old_letter:
        new_letter = random.choice(AMINOACIDS)
      return new_letter
    
    return ''.join(map(switch, sequence))
