from .Mutation import Mutation
import evodesign.Sequence as Sequence
import numpy as np





class Swap(Mutation):
  
  def _params(self) -> dict:
    params = super()._params()
    params['num_swaps'] = self._num_swaps
    return params
  


  def __init__(self, 
               mutation_prob: float = 1.0,
               num_swaps: int = 1
               ) -> None:
    """
    Mutation operation in which a given number of residues are randomly 
    chosen from the given amino acid sequence (with uniform probability and
    without replacement) and exchanged for a different amino acid.
    The replacing amino acid is chosen randomly with uniform probability. 

    For example, given the sequence 'AAAAAA', and assuming that only one 
    residue--the third one--is to be exchanged, then this operation would 
    produce the sequence'AADAAA'.

    Parameters
    ----------
    mutation_prob : float, optional
        The probability for applying this mutation operator over a given
        sequence. The default is 1.0 (i.e., mutate every and all sequences
        in a given population).
    num_swaps : int, optional
        The number of amino acids that will be exchanged in the sequence.
        The default is 1.
    """
    super().__init__(mutation_prob)
    self._num_swaps = num_swaps
  


  def mutate_sequence(self, 
                      sequence: str,
                      rng: np.random.Generator
                      ) -> str:
    """
    Modifies some amino acids residues in the given sequence.

    Parameters
    ----------
    sequence : str
        The amino acid sequence that will be modified. Each residue in this
        sequence must be represented by a single letter corresponding to one 
        of the 20 essential amino acids.
    rng : numpy.random.Generator
        The pseudo-random number generator.

    Returns
    -------
    str
        The modified sequence.
    """
    indices = rng.choice(len(sequence), self._num_swaps, replace=False)
    # in Python, strings are immutable
    seq_list = list(sequence)
    for i in indices:
      seq_list[i] = Sequence.swap_letter(rng, seq_list[i])
    return ''.join(seq_list)
