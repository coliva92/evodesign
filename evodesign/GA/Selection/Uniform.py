from .Selection import Selection
from ...Context import Context
import pandas as pd





class Uniform(Selection):
  
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
                                   replace=False)
    return population.loc[selection]
