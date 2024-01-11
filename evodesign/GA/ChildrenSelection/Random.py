from ChildrenSelection import ChildrenSelection
from ... import Population
import Choice





class Random(ChildrenSelection):

  @classmethod
  def name(cls) -> str:
    return 'GA_ChildrenSelection_Random'
  


  def __init__(self, 
               maxInputSize: int,
               bias: float = 0.5
               ) -> None:
    super().__init__(maxInputSize)
    self._bias = bias
    self._weights = ( bias, 1.0 - bias )
  


  def params_as_dict(self) -> dict:
    params = super().params_as_dict()
    params['bias'] = self._bias
    return params
  


  def select_children(self, children: Population) -> Population:
    # suponemos que la recombinación produjo dos hijos por cada par de padres
    children.individuals = [
      sister if Choice.flip_coin(self._weights) else brother
      for sister, brother in zip(children[0::2], children[1::2])
    ]
    return children
