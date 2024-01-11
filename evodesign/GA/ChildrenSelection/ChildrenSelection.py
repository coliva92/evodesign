from abc import ABC, abstractmethod, abstractclassmethod
from ... import Population





class ChildrenSelection(ABC):

  @abstractclassmethod
  def name(cls) -> str:
    raise NotImplementedError
  


  def __init__(self, maxInputSize: int) -> None:
    super().__init__()
    self._max_input_size = maxInputSize
  

  
  def params_as_dict(self) -> dict:
    return {}
  


  def __call__(self, children: Population) -> Population:
    # suponemos que en este punto ya se calculó la aptitud de los individuos, 
    # pero aun no están ordenados por aptitud
    children.individuals = children.individuals \
                           if len(children) <= self._max_input_size \
                           else children[:self._max_input_size]
    return self.select_children(children)
  


  @abstractmethod
  def select_children(self, children: Population) -> Population:
    raise NotImplementedError
