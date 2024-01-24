from abc import ABC, abstractclassmethod





class SettingsRetrievable(ABC):

  @abstractclassmethod
  def _name(cls) -> str:
    raise NotImplementedError
  


  def _params(self) -> dict:
    return {}
  


  def settings(self) -> dict:
    return { f'{self._name()}': self._params() }
