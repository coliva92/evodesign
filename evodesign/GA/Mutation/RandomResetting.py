from . import Mutation
import Sequence
from evodesign import Random as Choice





class RandomResetting(Mutation):

  @classmethod
  def name(cls) -> str:
    return 'GA_Mutation_RandomResetting'
  
  

  def __init__(self, 
               probability: float = 1.0, 
               residueProbability: float = 0.1
               ) -> None:
    super().__init__(probability)
    self._residue_prob = residueProbability
    self._weights = ( residueProbability, 1.0 - residueProbability )
  


  def params_as_dict(self) -> dict:
    params = super().params_as_dict()
    params['residueProbability'] = self._residue_prob
    return params



  def mutate_sequence(self, sequence: str) -> str:
    def switch(old_letter: str) -> str:
      if not Choice.coin_toss(self._weights):
        return old_letter
      return Sequence.switch_residue(old_letter)
    
    return ''.join(map(switch, sequence))
