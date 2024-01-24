from abc import ABC, abstractmethod
from ...SettingsRetrievable import SettingsRetrievable
import pandas as pd





class Selection(SettingsRetrievable, ABC):

  def _params(self) -> dict:
    return {
      'selectionSize': self._selection_size
    }
  

  
  def __init__(self, selectionSize: int) -> None:
    super().__init__()
    self._selection_size = selectionSize



  @abstractmethod
  def select_parents(self, population: pd.DataFrame) -> pd.DataFrame:
    raise NotImplementedError
  


  def __call__(self, population: pd.DataFrame) -> pd.DataFrame:
    """
    Selects a subset of individuals from the given population. Only those
    rows with a `True` value in the 'Survivor' column are considered.

    Parameters
    ----------
    population : pd.DataFrame
        The population to be sampled.

    Returns
    -------
    pd.DataFrame
        The selected subset of individuals.
    """
    survivors = population[population['Survivor']]
    return self.select_parents(survivors)
