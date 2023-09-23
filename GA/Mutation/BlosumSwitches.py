from .Mutation import Mutation
import random
import evodesign.Sequence as Sequence
import blosum as bl





class BlosumSwitches(Mutation):

  _matrix = bl.BLOSUM(62)
  _shifts = None



  @classmethod
  def name(cls) -> str:
    return 'GA_Mutation_BlosumSwitches'
  
  

  def __init__(self, 
               probability: float = 1.0, 
               numSwitches: int = 1
               ) -> None:
    super().__init__(probability)
    self._num_switches = numSwitches
    if not __class__._shifts:
      __class__._shifts = {}
      for residue in Sequence.AMINOACIDS:
        negatives = filter(lambda x: x[1] < 0 and x[0] not in { 'B', 'J', 'O', 'U', 'X', 'Z' }, 
                           [ (k, v) for k, v in __class__._matrix[residue].items() ])
        __class__._shifts[residue] = [ x[0] for x in negatives ]
  


  def params_json(self) -> dict:
    params = super().params_json()
    params['numSwitches'] = self._num_switches
    return params



  def mutate_sequence(self, sequence: str) -> str:
    positions = random.sample(range(len(sequence)), self._num_switches)
    # en Python, las cadenas son inmutables; hay que convertirlas a una lista
    seq_list = list(sequence)
    for i in positions:
      seq_list[i] = random.choice(__class__._shifts[seq_list[i]])
    return ''.join(seq_list)
