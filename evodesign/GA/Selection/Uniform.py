from .Selection import Selection
from ...Random import Random
import pandas as pd





class Uniform(Selection):

  @classmethod
  def _name(cls) -> str:
    return 'GA_Selection_Uniform'
  


  def __init__(self, selectionSize: int) -> None:
    """
    Random uniform sampling from the population without replacement.

    Parameters
    ----------
    selectionSize : int
        The number of individuals to be selected from the population.
    """
    super().__init__(selectionSize)
  


  def _select_parents(self, population: pd.DataFrame) -> pd.DataFrame:
    rng = Random.generator()
    selection = rng.choice(population.index, 
                           self._selection_size, 
                           replace=False)
    return population.iloc[selection]
