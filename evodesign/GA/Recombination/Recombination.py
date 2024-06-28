from abc import ABC, abstractmethod
from ...SettingsRetrievable import SettingsRetrievable
from typing import List
import evodesign.Population as Population
import evodesign.Random as Random
import pandas as pd





class Recombination(SettingsRetrievable, ABC):

  @abstractmethod
  def offspring_sequences(self, 
                          mother: str,
                          father: str
                          ) -> List[str]:
    raise NotImplementedError
  


  def _params(self) -> dict:
    return { 'probability': self._probability }
  


  def __init__(self, probability: float) -> None:
    super().__init__()
    self._probability = probability
  


  def __call__(self, 
               parents: pd.DataFrame,
               generationId: int = 0
               ) -> pd.DataFrame:
    """
    Produces a table of new sequences generated by mixing the sequences 
    of even rows with those of the odd rows.

    Parameters
    ----------
    parents : pandas.DataFrame
        The table which rows will be mixed to produce the new table. 
        It is assumed that there is an even number of rows
        in this table, otherwise the last row is dropped.
    generationId : int, optional
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
      children += self.offspring_sequences(mother, father) \
                  if Random.coin_toss(self._probability) \
                  else [ mother, father ]
    return Population.create(children, generationId)
  