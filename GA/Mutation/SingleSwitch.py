from .Mutation import Mutation
import random
import evodesign.Sequence as Sequence





class SingleSwitch(Mutation):

  @classmethod
  def get_name(cls) -> str:
    return 'GA_Mutation_SingleSwitch'



  def mutate_sequence(self, sequence: str) -> str:
    i = random.randint(0, len(sequence) - 1)
    # en Python, las cadenas son inmutables; hay que convertirlas a una lista 
    seq_list = list(sequence)
    seq_list[i] = Sequence.switch_residue(seq_list[i])
    return ''.join(seq_list)
