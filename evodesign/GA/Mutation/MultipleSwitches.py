from . import Mutation
import random
import Sequence





class MultipleSwitches(Mutation):

  @classmethod
  def name(cls) -> str:
    return 'GA_Mutation_MultipleSwitches'
  
  

  def __init__(self, 
               probability: float = 1.0, 
               numSwitches: int = 1
               ) -> None:
    super().__init__(probability)
    self._num_switches = numSwitches
  


  def params_as_dict(self) -> dict:
    params = super().params_as_dict()
    params['numSwitches'] = self._num_switches
    return params



  def mutate_sequence(self, sequence: str) -> str:
    positions = random.sample(range(len(sequence)), self._num_switches)
    # en Python, las cadenas son inmutables; hay que convertirlas a una lista
    seq_list = list(sequence)
    for i in positions:
      seq_list[i] = Sequence.switch_residue(seq_list[i])
    return ''.join(seq_list)
