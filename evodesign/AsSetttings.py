from abc import ABC, abstractclassmethod





class AsSettings(ABC):

  @abstractclassmethod
  def _name(cls) -> str:
    raise NotImplementedError
  


  def _params(self) -> dict:
    return {}
  


  def as_settings(self) -> dict:
    return { f'{self._name()}': self._params() }
