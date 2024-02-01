from .Recombination import Recombination
from typing import List
from ...Random import Random





class SinglePointCrossover(Recombination):

  @classmethod
  def _class_name(cls) -> str:
    return 'GA.Recombination.SinglePointCrossover'



  def __init__(self) -> None:
    """
    Creates new sequences by splitting the parent sequences at a randomly 
    chosen residue position and mixing one part from one of the parents with
    the complementary part from the other.

    For example, given the sequences 'AAAAAA' and 'DDDDDD', and assuming that 
    it was chosen to split these sequences at the 4th position, then this 
    operator would produce the sequences 'AAAADD' and 'DDDDAA'.
    """
    super().__init__()
  


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
    i = rng.integers(0, len(mother))
    sister = mother[0:i] + father[i:]
    brother = father[0:i] + mother[i:]
    return [ sister, brother ]
