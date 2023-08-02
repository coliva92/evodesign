from .Mutation import Mutation
import random
from evodesign.Sequence import AMINOACIDS





class SingleSwitch(Mutation):

  @classmethod
  def get_name(cls) -> str:
    return 'GA_Mutation_SingleSwitch'



  def mutate_sequence(self, sequence: str) -> str:
    i = random.randint(0, len(sequence) - 1)
    seq_list = list(sequence)
    original_letter = seq_list[i]
    while sequence[i] == original_letter:
      seq_list[i] = random.choice(AMINOACIDS)
    return ''.join(seq_list)
