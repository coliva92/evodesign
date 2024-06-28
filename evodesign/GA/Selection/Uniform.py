from .Selection import Selection
import evodesign.Random as Random
import pandas as pd





class Uniform(Selection):
  
  def select_parents(self, population: pd.DataFrame) -> pd.DataFrame:
    """
    Selects a subset of individuals from the given population.

    Parameters
    ----------
    population : pandas.DataFrame
        The population to be sampled.

    Returns
    -------
    pandas.DataFrame
        The selected subset of individuals.
    """
    rng = Random.generator()
    selection = rng.choice(population.index, 
                           len(population), 
                           replace=False)
    return population.loc[selection]
