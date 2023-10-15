from .ChildrenSelection import ChildrenSelection
from evodesign import Population
import evodesign.Choice as Choice





class BetterFitness(ChildrenSelection):

  @classmethod
  def name(cls) -> str:
    return 'GA_ChildrenSelection_BetterFitness'
  


  def __init__(self, maxInputSize: int, probability: float = 1.0) -> None:
    super().__init__(maxInputSize)
    self._probability = probability
    self._weights = ( probability, 1.0 - probability )
  


  def params_json(self) -> dict:
    params = super().params_json()
    params['probability'] = self._probability
    return params
  

  
  def select_children(self, children: Population) -> Population:
    # suponemos que la recombinación produjo dos hijos por cada par de padres
    temp_children = []
    for sister, brother in zip(children[0::2], children[1::2]):
      if Choice.flip_coin(self._weights):
        temp_children.append(sister if sister >= brother else brother)
      else:
        temp_children.append(sister if sister < brother else brother)
    children.individuals = temp_children
    return children
