from .ChildSelection import ChildSelection
from evodesign import Population





class BestChild(ChildSelection):

  @classmethod
  def name(cls) -> str:
    return 'GA_ChildSelection_SingleBest'
  

  
  def select_children(self, children: Population) -> Population:
    # suponemos que la recombinación produjo dos hijos por cada par de padres
    children.individuals = [
      sister if sister > brother else brother
      for sister, brother in zip(children[0::2], children[1::2])
    ]
    return children
