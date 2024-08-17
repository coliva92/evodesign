from abc import ABC, abstractmethod
from ...SettingsRetrievable import SettingsRetrievable
from ...Context import Context
import evodesign.Utils as Utils
import pandas as pd
from typing import List





class Selection(SettingsRetrievable, ABC):
  
    def _params(self) -> dict:
        params = super()._params()
        params["two_offspring"] = self._two_offspring
  


    def __init__(self, two_offspring: bool = True) -> None:
        super().__init__()
        self._two_offspring = two_offspring
    


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
        num_couple_selections = len(survivors) // 2 \
                                if self._two_offspring \
                                else len(survivors)
        parents = pd.DataFrame(columns=survivors.columns)
        for _ in range(num_couple_selections):
            selection = self.select_parent_couple(population, context)
            parents = Utils.df_append(parents, selection[0])
            parents = Utils.df_append(parents, selection[1])
        return parents
  


    @abstractmethod
    def select_parent_couple(self, 
                             population: pd.DataFrame,
                             context: Context
                             ) -> List[pd.Series]:
        raise NotImplementedError
