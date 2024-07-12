from abc import ABC, abstractmethod
from ...SettingsRetrievable import SettingsRetrievable
import evodesign.Utils as Utils
import numpy as np
import pandas as pd
from typing import Optional





class Mutation(SettingsRetrievable, ABC):

  def __init__(self, mutation_prob: float = 1.0) -> None:
    super().__init__()
    self._mutation_prob = mutation_prob



  def _params(self) -> dict:
    return {
      'mutation_prob': self._mutation_prob
    }



  def __call__(self, 
               rng: np.random.Generator,
               children: pd.DataFrame
               ) -> pd.DataFrame:
    """
    Modifies the amino acid sequences of a select subset of rows in the given 
    table.

    Parameters
    ----------
    rng : numpy.random.Generator
        The pseudo-random number generator.
    children : pandas.DataFrame
        The table from which a subset will be selected and modified.
    """
    indices = children.apply(lambda row: Utils.coin_toss(rng, self._mutation_prob), 
                             axis=1)
    children.loc[indices, 'sequence'] = \
      children.loc[indices, 'sequence'].apply(self.mutate_sequence, args=(rng))
    return children
  


  @abstractmethod
  def mutate_sequence(self, 
                      sequence: str,
                      rng: np.random.Generator
                      ) -> str:
    raise NotImplementedError
