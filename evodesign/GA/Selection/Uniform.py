from Selection import Selection
from ...Random import Random
import pandas as pd





class Uniform(Selection):

  @classmethod
  def _class_name(cls) -> str:
    return 'GA.Selection.Uniform'
  


  def __init__(self, numSelectedCouples: int) -> None:
    """
    Random uniform sampling from the population without replacement.

    Parameters
    ----------
    numSelectedCouples : int
        The number of parent couples to be selected from the population.
    """
    super().__init__(numSelectedCouples)
  


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
                           self._selection_size, 
                           replace=False)
    return population.loc[selection]
