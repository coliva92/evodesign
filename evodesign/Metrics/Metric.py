from abc import ABC, abstractmethod
from ..SettingsRetrievable import SettingsRetrievable
from typing import Dict





class Metric(SettingsRetrievable, ABC):

  @abstractmethod
  def column_name(self) -> str:
    raise NotImplementedError
  


  @abstractmethod
  def compute_value(self, **kwargs) -> float:
    raise NotImplementedError
  


  def __call__(self, **kwargs) -> float:
    """
    Computes the specified metric for a given model backbone or sequence.

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
    otherMetrics : Dict[str, int  |  float  |  str], optional
        The previously computed values for other metrics. Default is an empty 
        dictionary.

    Returns
    -------
    float
        The final metric value.
    """
    if 'otherMetrics' not in kwargs:
      kwargs['otherMetrics'] = {}
    if self.column_name() in kwargs['otherMetrics']:
      return kwargs['otherMetrics'][self.column_name()]
    return self.compute_value(**kwargs)
