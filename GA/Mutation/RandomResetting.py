from .Mutation import Mutation
import random
from evodesign.Sequence import AMINOACIDS





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
    def mutate_residue(residue: str) -> str:
      result = residue
      if random.choices(Mutation._options, self._weights)[0]:
        while result == residue:
          result = random.choice(AMINOACIDS)
      return result
    return ''.join(map(mutate_residue, sequence))
