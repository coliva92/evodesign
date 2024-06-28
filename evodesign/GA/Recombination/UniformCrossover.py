from .Recombination import Recombination
from typing import List
import evodesign.Random as r





class UniformCrossover(Recombination):

  @classmethod
  def _class_name(cls) -> str:
    return 'GA.Recombination.UniformCrossover'
  


  def _params(self) -> dict:
    params = super()._params()
    params['maskBias'] = self._parent_bias
    return params
  
  

  def __init__(self,
               maskBias: float = 0.5,
               probability: float = 1.0
               ) -> None:
    """
    Randomly generates a binary mask of the same length as the parent sequences
    and uses it to construct a new sequence, residue by residue, by choosing
    one amino acid from one parent or the other for the corresponding position
    in the generated mask.

    For example, given the sequences 'AAAAAA' and 'DDDDDD', and assuming the 
    random binary mask was 011010, then this operator would produce the 
    sequences 'ADDADA' and 'DAADAD'.

    Parameters
    ----------
    maskBias : float, optional
        The probability for generating one binary value over the other in the
        random mask. The default is 0.5.
    """
    super().__init__(probability)
    self._weights = ( maskBias, 1.0 - maskBias )
    self._parent_bias = maskBias



  def offspring_sequences(self, 
                          mother: str,
                          father: str
                          ) -> List[str]:
    """
    Mix partial information from the two given sequences to create two new
    sequences. It is assumed that both sequences have equal length.

    In both sequences, each residue must be represented by a single letter
    corresponding to one of the 20 essential amino acids.

    Parameters
    ----------
    mother : str
        One of the sequences to be mixed.
    father : str
        The other sequence to be mixed.

    Returns
    -------
    List[str]
        The two sequences produced.
    """
    rng = r.generator()
    mask = rng.choice([ 0, 1 ], size=len(mother), p=self._weights)
    parents = ( mother, father )
    sister = ''.join([ parents[p][i] for i, p in enumerate(mask) ])
    parents = ( father, mother )
    brother = ''.join([ parents[p][i] for i, p in enumerate(mask) ])
    return [ sister, brother ]
