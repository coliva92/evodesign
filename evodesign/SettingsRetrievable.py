from abc import ABC, abstractmethod





class SettingsRetrievable(ABC):

  @classmethod
  @abstractmethod
  def _class_name(cls) -> str:
    raise NotImplementedError
  


  def _params(self) -> dict:
    return {}
  


  def settings(self) -> dict:
    return { f'{self._class_name()}': self._params() }
