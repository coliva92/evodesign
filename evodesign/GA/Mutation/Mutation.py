from abc import ABC, abstractmethod
from ...SettingsRetrievable import SettingsRetrievable
import evodesign.Random as Random
import pandas as pd





class Mutation(SettingsRetrievable, ABC):

  def __init__(self, mutProb: float = 1.0) -> None:
    super().__init__()
    self._mutation_prob = mutProb



  def _params(self) -> dict:
    return {
      'mutProb': self._mutation_prob
    }


  
  @abstractmethod
  def mutate_sequence(self, sequence: str) -> str:
    raise NotImplementedError



  def __call__(self, children: pd.DataFrame) -> pd.DataFrame:
    """
    Modifies the amino acid sequences of a select subset of rows in the given 
    table.

    Parameters
    ----------
    children : pandas.DataFrame
        The table from which a subset will be selected and modified.
    """
    indices = children.apply(lambda row: Random.coin_toss(self._mutation_prob), 
                             axis=1)
    children.loc[indices, 'sequence'] = \
      children.loc[indices, 'sequence'].apply(self.mutate_sequence)
    return children
