from .Mutation import Mutation
import evodesign.Sequence as Sequence
import evodesign.Utils as Utils
import numpy as np





class RandomResetting(Mutation):
  
  def _params(self) -> dict:
    params = super()._params()
    params['swap_prob'] = self._swap_prob
    return params
  
  

  def __init__(self, 
               mutation_prob: float = 1.0,
               swap_prob: float = 0.1
               ) -> None:
    """
    Mutation operation in which a binary mask of the same length as the given 
    sequence is randomly produced, and where all the residues in which the 
    corresponding position in said mask is 1 are exchanged for a different 
    amino acid.

    For example, given the sequence 'AAAAAA', and assuming the mask 110001 was
    randomly generated, then this operation would produce the sequence
    'DFAAAL'.

    Parameters
    ----------
    mutation_prob : float, optional
        The probability for applying this mutation operator over a given
        sequence. The default is 1.0.
    swap_prob : float, optional
        The probability for each residue in the sequence to be exchanged
        for another amino acid. The default is 0.1.
    """
    super().__init__(mutation_prob)
    self._swap_prob = swap_prob
    self._residue_weights = ( swap_prob, 1.0 - swap_prob )
  


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
    return ''.join([
      Sequence.swap_letter(rng, x) if Utils.coin_toss(rng, self._swap_prob) else x
      for x in sequence
    ])
