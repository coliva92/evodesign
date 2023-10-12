from .ChildSelection import ChildSelection
from evodesign import Population
import evodesign.Choice as Choice





class RandomChild(ChildSelection):

  @classmethod
  def name(cls) -> str:
    return 'GA_ChildSelection_SingleRandom'
  


  def __init__(self, 
               maxInputSize: int,
               sisterProbability: float = 0.5
               ) -> None:
    super().__init__(maxInputSize)
    self._sister_probability = sisterProbability
    self._weights = ( sisterProbability, 1.0 - sisterProbability )
  


  def params_json(self) -> dict:
    params = super().params_json()
    params['sisterProbability'] = self._sister_probability
  


  def select_children(self, children: Population) -> Population:
    # suponemos que la recombinación produjo dos hijos por cada par de padres
    children.individuals = [
      sister if Choice.flip_coin(self._weights) else brother
      for sister, brother in zip(children[0::2], children[1::2])
    ]
    return children
