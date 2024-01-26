from Mutation import Mutation
from ...Random import Random
from ...Sequence import Sequence





class Switch(Mutation):

  @classmethod
  def _class_name(cls) -> str:
    return 'GA.Mutation.Switch'
  


  def _params(self) -> dict:
    params = super()._params()
    params['numSwitches'] = self._num_switches
    return params
  


  def __init__(self, 
               mutProb: float = 1.,
               numSwitches: int = 1
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
    mutProb : float, optional
        The probability for applying this mutation operator over a given
        sequence. The default is 1.0 (i.e., mutate every and all sequences
        in a given population).
    numSwitches : int, optional
        The number of amino acids that will be exchanged in the sequence.
        The default is 1.
    """
    super().__init__(mutProb)
    self._num_switches = numSwitches
  


  def mutate_sequence(self, sequence: str) -> str:
    """
    Modifies some amino acids residues in the given sequence.

    Parameters
    ----------
    sequence : str
        The amino acid sequence that will be modified. Each residue in this
        sequence must be represented by a single letter corresponding to one 
        of the 20 essential amino acids.

    Returns
    -------
    str
        The modified sequence.
    """
    rng = Random.generator()
    indices = rng.choice(len(sequence), self._num_switches, replace=False)
    # in Python, strings are immutable
    seq_list = list(sequence)
    for i in indices:
      seq_list[i] = Sequence.switch_residue(seq_list[i])
    return ''.join(seq_list)
