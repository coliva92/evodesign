from .Selection import Selection
import evodesign.Random as r
import pandas as pd





class Uniform(Selection):

  @classmethod
  def _class_name(cls) -> str:
    return 'GA.Selection.Uniform'
  


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
    rng = r.generator()
    selection = rng.choice(population.index, 
                           len(population), 
                           replace=False)
    return population.loc[selection]
