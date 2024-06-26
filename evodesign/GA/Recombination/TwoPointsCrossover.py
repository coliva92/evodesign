from .Recombination import Recombination
from typing import List
import evodesign.Random as Random





class TwoPointsCrossover(Recombination):

  def __init__(self, probability: float = 1.0) -> None:
    """
    Creates new sequences by splitting the parent sequences at three randomly
    chosen residue positions and mixing two parts from one of the parents with
    one complementary part from the other.

    For example, given the sequences 'AAAAAA' and 'DDDDDD', and assuming that
    these sequences are to be split at the 2nd and 2nd to last positions, then
    this operator would produce the sequences 'AADDDA' and 'DDAAAD'.
    """
    super().__init__(probability)
  


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
    rng = Random.generator()
    n = len(mother)
    i = rng.integers(0, n) 
    j = rng.integers(i, n)
    sister = mother[0:i] + father[i:j] + mother[j:]
    brother = father[0:i] + mother[i:j] + father[j:]
    return [ sister, brother ]
