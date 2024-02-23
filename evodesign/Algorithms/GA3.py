from typing import List
from Bio.PDB.Atom import Atom
from .GenericGA import GenericGA
from ..Statistics import Statistics
import evodesign.Operators as op
import pandas as pd





class GA3(GenericGA):

  @classmethod
  def _class_name(cls) -> str:
    return 'Algorithms.GA3'
  


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
                  children: pd.DataFrame, 
                  reference: List[Atom]
                  ) -> pd.DataFrame:
    elite = population.iloc[0].copy()
    elite['survivor'] = False
    elite['generation_id'] = children.iloc[0]['generation_id']
    mixed = pd.concat([ elite.to_frame().T, children ], ignore_index=True)
    for i, row in mixed.iloc[1:].iterrows():
      if op.comes_first(elite, row, self._sort_cols, self._sort_ascending):
        continue
      # put winning row in the elite position
      mixed.iloc[0] = row
      # move last row to the winning row's original position
      mixed.iloc[i] = mixed.iloc[-1]
      # put elite solution (which is no longer the elite) in the last position
      mixed.iloc[-1] = elite
      break
    # if last pop's elite solution wins, it is preserved;
    # if it loses to another solution in children, it is removed and all
    # children survive to the next generation
    mixed.loc[:self._pop_size - 1, 'survivor'] = True
    return mixed
