from abc import ABC, abstractmethod
from ...SettingsRetrievable import SettingsRetrievable
import pandas as pd
from ...Context import Context





class Selection(SettingsRetrievable, ABC):
  
  def __call__(self,
               population: pd.DataFrame,
               context: Context
               ) -> pd.DataFrame:
    """
    Selects a subset of individuals from the given population. Only those
    rows with a `True` value in the 'survivor' column are considered.

    Parameters
    ----------
    population : pandas.DataFrame
        The population to be sampled.
    context : Context
        The context data used by the calling evolutionary algorithm.

    Returns
    -------
    pandas.DataFrame
        The selected subset of individuals.
    """
    survivors = population[population['survivor']]
    return self.select_parents(survivors, context)
  


  @abstractmethod
  def select_parents(self, 
                     population: pd.DataFrame,
                     context: Context
                     ) -> pd.DataFrame:
    raise NotImplementedError
