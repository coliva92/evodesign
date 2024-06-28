from .GenericGA import GenericGA
from ..Statistics import Statistics
import evodesign.Utils as Utils
import pandas as pd





class GA1(GenericGA):

  @classmethod
  def _class_name(cls) -> str:
    return 'Algorithms.GA1'
  


  def compute_statistics(self, population: pd.DataFrame) -> pd.Series:
    top_solution = population.iloc[0].copy()
    survivors = population[population['survivor']]
    top_solution['sequence_identity'] = \
      Statistics.average_sequence_identity(survivors)
    top_solution['lost_amino_acids'] = \
      Statistics.average_amino_acid_loss(survivors)
    return top_solution
  


  def replacement(self, 
                  population: pd.DataFrame, 
                  children: pd.DataFrame
                  ) -> pd.DataFrame:
    children.sort_values(by=self._sort_cols, 
                         ascending=self._sort_ascending, 
                         inplace=True,
                         ignore_index=True)
    last_pop = population[population['survivor']].copy()
    last_pop['survivor'] = False
    last_pop['generation_id'] = children.iloc[0]['generation_id']
    mixed = Utils.merge(last_pop, 
                        children, 
                        self._sort_cols, 
                        self._sort_ascending)
    mixed.loc[:self._pop_size - 1, 'survivor'] = True
    return mixed
