from .Mutation import Mutation
import random
from evodesign.Sequence import AMINOACIDS





class MultipleSwitches(Mutation):

  @classmethod
  def get_name(cls) -> str:
    return 'GA_Mutation_MultipleSwitches'
  
  

  def __init__(self, 
               probability: float = 1.0, 
               numSwitches: int = 1
               ) -> None:
    super().__init__(probability)
    self._num_switches = numSwitches
  


  def as_dict(self) -> dict:
    params = super().as_dict()
    params['numSwitches'] = self._num_switches
    return params



  def mutate_sequence(self, sequence: str) -> str:
    positions = random.sample(list(range(len(sequence))), self._num_switches)
    # en Python, las cadenas son inmutables; hay que convertirlas a una lista
    seq_list = list(sequence)
    for i, original_letter in zip(positions, seq_list):
      while seq_list[i] == original_letter: 
        seq_list[i] = random.choice(AMINOACIDS)
    return ''.join(seq_list)
