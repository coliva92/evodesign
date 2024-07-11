from .Metric import Metric
from typing import Optional





class Cyclization(Metric):

  def __init__(self, column: Optional[str] = None) -> None:
    super().__init__(column)



  def compute_value(self, **kwargs) -> float:
    # distance between the N atom of the first residue and the C 
    # atom of the last
    return kwargs['model'][0] - kwargs['model'][-2]
