from abc import ABC, abstractmethod
from typing import List
from ..Metrics.Metric import Metric
from ..SettingsRetrievable import SettingsRetrievable
import pandas as pd





class FitnessFunction(SettingsRetrievable, ABC):

  @classmethod
  @abstractmethod
  def column_name(cls) -> str:
    raise NotImplementedError
  


  def _params(self) -> dict:
    return { 'upperBound': self._upper_bound }
  


  def __init__(self, 
               upperBound: float, 
               metrics: List[Metric]
               ) -> None:
    super().__init__()
    self._upper_bound = upperBound
    self._metrics = metrics
  


  def upper_bound(self) -> float:
    """
    Returns
    -------
    float
        The highest value the fitness function can achieve before triggering
        the termination of the evolutionary algorithm.
    """
    return self._upper_bound


  
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
        The backbone of the target protein. The fitness of each individual
        in the population can be computed from this backbone.
    refSequence: str, optional
        The amino acid sequence of the target protein. The fitness of each
        individual in the population can be computed from this sequence.
        Default is `None`.
    sequence : str, optional
        The amino acid sequence for which the fitness will be computed.
        Each residue must be represented by a single letter corresponding to
        one of the 20 essential amino acids.
    sequence_id : str, optional
        The unique identifier for the given sequence.
    plddt : float, optional
        The predicted lDDT value computed by the protein structure prediction
        algorithm that was used for obtaining the given model backbone.

    Returns
    -------
    pandas.Series
        The computed metrics and fitness value for the given model or sequence.
    """
    values = {}
    for metric in self._metrics:
      values = { **values, **metric(**kwargs) }
    values[self.column_name()] = self.compute_fitness(**{ **kwargs, **values })
    return pd.Series(values)
