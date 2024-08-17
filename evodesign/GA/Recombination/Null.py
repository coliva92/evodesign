from .Recombination import Recombination
from typing import List
import evodesign.Population as Population
import evodesign.Utils as Utils
import pandas as pd
import numpy as np





class Null(Recombination):

  def __call__(self, 
               rng: np.random.Generator,
               parents: pd.DataFrame,
               generation_id: int = 0
               ) -> pd.DataFrame:
    children = parents.copy()
    children["generation_id"] = generation_id
    return children
  


  def offspring_sequences(self, 
                          rng: np.random.Generator,
                          mother: str,
                          father: str
                          ) -> List[str]:
    raise NotImplementedError
