from .Selection import Selection
from ...Context import Context
import pandas as pd





class Uniform(Selection):
  
    def _params(self) -> dict:
       params = super()._params()
       params["replacement"] = self._replacement
       return params
    


    def __init__(self, replacement: bool = False) -> None:
        """
        Selection operation in which a subset of population individuals is choosed
        randomly with uniform probability. 

        Parameters
        ----------
        replace : bool, optional
            Indicates if the random selection should be performed with or without 
            replacement. Default is `False`.
        """
        super().__init__()
        self._replacement = replacement
    


    def select_parents(self, 
                       population: pd.DataFrame,
                       context: Context
                       ) -> pd.DataFrame:
      """
      Selects a subset of individuals from the given population.
  
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
      selection = context.rng.choice(population.index, 
                                     len(population), 
                                     replace=self._replacement)
      return population.loc[selection]
