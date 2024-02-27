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
    children.sort_values(by=self._sort_cols, 
                         ascending=self._sort_ascending,
                         inplace=True,
                         ignore_index=True)
    if op.comes_first(elite, children.iloc[0], self._sort_cols, self._sort_ascending):
      mixed = pd.concat([ elite.to_frame().T, children ], ignore_index=True)
      mixed['survivor'] = True
      mixed.loc[len(mixed) - 1, 'survivor'] = False
    elif op.comes_first(children.iloc[0], elite, self._sort_cols, self._sort_ascending):
      children['survivor'] = True
      mixed = pd.concat([ children, elite.to_frame().T ], ignore_index=True)
    return mixed
