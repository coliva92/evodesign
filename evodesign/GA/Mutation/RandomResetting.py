from Mutation import Mutation
from ...Sequence import Sequence
from ...Random import Random





class RandomResetting(Mutation):

  @classmethod
  def _name(cls) -> str:
    return 'GA.Mutation.RandomResetting'
  


  def _params(self) -> dict:
    params = super()._params()
    params['exchangeProb'] = self._exchange_prob
    return params
  
  

  def __init__(self, 
               mutProb: float = 1.0, 
               exchangeProb: float = 0.1
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
    mutProb : float, optional
        The probability for applying this mutation operator over a given
        sequence. The default is 1.0 (i.e., mutate every and all sequences
        in a given population).
    exchangeProb : float, optional
        The probability for each residue in the sequence to be exchanged
        for another amino acid. 
        The default is 0.1.
    """
    super().__init__(mutProb)
    self._exchange_prob = exchangeProb
    self._residue_weights = ( exchangeProb, 1. - exchangeProb )
  


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
    return ''.join([
      Sequence.switch_residue(x) if Random.coin_toss(self._exchange_prob) else x
      for x in sequence
    ])
