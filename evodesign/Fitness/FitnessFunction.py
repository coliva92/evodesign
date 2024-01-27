from abc import ABC, abstractmethod, abstractclassmethod
from typing import Dict, List, Optional
from Bio.PDB.Atom import Atom
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
  def compute_fitness(self, 
                      metricValues: Dict[str, float],
                      sequence: Optional[str] = None
                      ) -> float:
    raise NotImplementedError



  def __call__(self,
               model: List[Atom], 
               reference: List[Atom],
               sequence: Optional[str] = None
               ) -> pd.Series:
    values = {
      metric.column_name(): metric(model, reference, sequence)
      for metric in self._metrics
    }
    values[self.column_name()] = self.compute_fitness(metricValues=values, 
                                                      sequence=sequence)
    return pd.Series(values)
