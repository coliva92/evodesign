from abc import ABC, abstractmethod, abstractclassmethod





class Metric(ABC):

  @abstractclassmethod
  def column_name(cls) -> str:
    raise NotImplemented
  

  
  @abstractmethod
  def __call__(self, **kwargs) -> float:
    raise NotImplementedError
