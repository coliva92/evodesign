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
    values = {
      metric.column_name(): metric(**kwargs)
      for metric in self._metrics
    }
    values[self.column_name()] = self.compute_fitness(**kwargs)
    return pd.Series(values)
