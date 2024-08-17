from .Selection import Selection
from ...Context import Context
import pandas as pd
from typing import List





class Null(Selection):

    def __call__(self,
                 population: pd.DataFrame,
                 context: Context
                 ) -> pd.DataFrame:
        survivors = population[population['survivor']]
        return survivors



    def select_parent_couple(self, 
                             population: pd.DataFrame,
                             context: Context
                             ) -> List[pd.Series]:
        raise NotImplementedError
