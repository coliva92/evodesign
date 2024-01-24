from abc import ABC, abstractmethod
from ...SettingsRetrievable import SettingsRetrievable
from ...Random import Random
import pandas as pd





class Mutation(SettingsRetrievable, ABC):

  def __init__(self, mutProb: float = 1.0) -> None:
    super().__init__()
    self._weights = ( mutProb, 1.0 - mutProb )
    self._mutation_prob = mutProb



  def _params(self) -> dict:
    return {
      'mutProb': self._mutation_prob
    }


  
  @abstractmethod
  def mutate_sequence(self, sequence: str) -> str:
    raise NotImplementedError



  def __call__(self, children: pd.DataFrame) -> None:
    """
    Modifies the amino acid sequences of a select subset of the given 
    collection of individuals.

    Parameters
    ----------
    children : pd.DataFrame
        The collection of individuals from which a subset will be selected 
        and modified.
    """
    indices = children.apply(lambda row: Random.coin_toss(self._weights), 
                             axis=1)
    children[indices, 'Sequence'].apply(self.mutate_sequence, inplace=True)
