from abc import ABC, abstractmethod, abstractclassmethod
from typing import List
from ..Metrics.Metric import Metric
from ..SettingsRetrievable import SettingsRetrievable
import pandas as pd





class FitnessFunction(SettingsRetrievable, ABC):

  @abstractclassmethod
  def column_name(cls) -> str:
    raise NotImplementedError



  @abstractclassmethod
  def upper_bound(cls) -> float:
    raise NotImplementedError
  


  def __init__(self, metrics: List[Metric]) -> None:
    super().__init__()
    self._metrics = metrics


  
  @abstractmethod
  def compute_fitness(self, **kwargs) -> float:
    raise NotImplementedError



  def __call__(self, **kwargs) -> pd.Series:
    """
    Computes the fitness for a given model backbone or sequence.

    Parameters
    ----------
    model : List[Bio.PDB.Atom.Atom]
        The model backbone for which the fitness will be computed.
    reference : List[Bio.PDB.Atom.Atom]
        The reference backbone. The model usually gets compared against the
        backbone in order to compute the fitness value of the former.
    sequence : str
        The name amino acid sequence for which the fitness will be computed.
        Each residue must be represented by a single letter corresponding to
        one of the 20 essential amino acids.

    Returns
    -------
    pandas.Series
        The computed metrics and fitness value for the given model or sequence.
    """
    values = {
      metric.column_name(): metric(**kwargs)
      for metric in self._metrics
    }
    values[self.column_name()] = self.compute_fitness(**{ **kwargs, **values })
    return pd.Series(values)
