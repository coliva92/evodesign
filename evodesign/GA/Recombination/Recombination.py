from abc import ABC, abstractmethod
from ...SettingsRetrievable import SettingsRetrievable
from typing import List
import evodesign.Population as Population
import evodesign.Utils as Utils
import pandas as pd
import numpy as np





class Recombination(SettingsRetrievable, ABC):

  def _params(self) -> dict:
    params = super()._params()
    params["two_offspring"] = self._two_offspring
    return params
  


  def __init__(self, two_offspring: bool = True) -> None:
    super().__init__()
    self._two_offspring = two_offspring
  

  
  def __call__(self, 
               rng: np.random.Generator,
               parents: pd.DataFrame,
               generation_id: int = 0
               ) -> pd.DataFrame:
    """
    Produces a table of new sequences generated by mixing the sequences 
    of even rows with those of the odd rows.

    Parameters
    ----------
    rng : numpy.random.Generator
        The pseudo-random number generator.
    parents : pandas.DataFrame
        The table which rows will be mixed to produce the new table. 
        It is assumed that there is an even number of rows
        in this table, otherwise the last row is dropped.
    generation_id : int, optional
        The unique identifier for the population being created. 
        The default value is 0.

    Returns
    -------
    pandas.DataFrame
        The table of the produced sequences.
    """
    if len(parents) % 2 != 0:
      parents = parents.iloc[:-1]
    children = []
    for i in range(0, len(parents), 2):
      mother = parents.iloc[i]['sequence']
      father = parents.iloc[i + 1]['sequence']
      temp = self.offspring_sequences(rng, mother, father)
      if self._two_offspring:
        children += temp
      else:
        children.append(temp[Utils.coin_toss(rng)])
    return Population.create(children, generation_id)
  


  @abstractmethod
  def offspring_sequences(self, 
                          rng: np.random.Generator,
                          mother: str,
                          father: str
                          ) -> List[str]:
    raise NotImplementedError
  