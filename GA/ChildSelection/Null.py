from .ChildSelection import ChildSelection
from evodesign import Population





class Null(ChildSelection):

  @classmethod
  def name(cls) -> str:
    return 'GA_ChildSelection_None'



  def select_children(self, children: Population) -> Population:
    return children
  