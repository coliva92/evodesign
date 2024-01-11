from ChildrenSelection import ChildrenSelection
from ... import Population





class Null(ChildrenSelection):

  @classmethod
  def name(cls) -> str:
    return 'GA_ChildrenSelection_None'



  def select_children(self, children: Population) -> Population:
    return children
  