from abc import ABC
from Utils import EvoDesignTestCase
from evodesign.GA.Selection.Overselection import Overselection
from evodesign.Population import Population
import pandas as pd





class SelectionTests(EvoDesignTestCase, ABC):

  def setUp(self):
    self.sequences = [
      'AAAAAAAA',
      'CCCCCCCC',
      'DDDDDDDD',
      'EEEEEEEE',
      'FFFFFFFF',
      'GGGGGGGG',
      'HHHHHHHH',
      'IIIIIIII'
    ]
    self.pop = Population.create(self.sequences)
    self.pop['survivor'] = True
    self.seq_length = len(self.sequences[0])
    self.num_couples = 3
    self.selection_size = 6
  


  def couples_are_different_sequences(self, parents):
    num_couples_diff_seq = sum([ 
      parents.iloc[i]['sequence'] != parents.iloc[i + 1]['sequence'] 
      for i in range(0, len(parents), 2) 
    ])
    self.assertEqual(num_couples_diff_seq, self.num_couples)
    return True



  def valid_couples(self, parents):
    self.assertIsInstance(parents, pd.DataFrame)
    parents_length = len(parents)
    self.assertEqual(parents_length, self.selection_size)
    self.couples_are_different_sequences(parents)





class OverselectionTests(SelectionTests):

  def setUp(self):
    super().setUp()
    self.upper_size = 3
    self.upper_bin = dict.fromkeys([ seq for seq in self.sequences[:self.upper_size] ])
    self.lower_bin = dict.fromkeys([ seq for seq in self.sequences[self.upper_size:] ])
  
  

  def sequences_from_bin(self, parents, correctBin):
    num_upper_in_upper_parents = sum([
      row['sequence'] in correctBin 
      for _, row in parents.iterrows() 
    ])
    self.assertEqual(num_upper_in_upper_parents, self.selection_size)



  def test_upper_bin_selection(self):
    selection = Overselection(numCouples=self.num_couples, 
                              upperSize=self.upper_size, 
                              upperProb=1.0, 
                              lowerProb=0.0)
    parents = selection(self.pop)
    self.valid_population(parents, 
                          self.selection_size, 
                          self.seq_length, 
                          survivorFlag=True)
    self.valid_couples(parents)
    self.sequences_from_bin(parents, self.upper_bin)
  


  def test_lower_bin_selection(self):
    selection = Overselection(numCouples=self.num_couples, 
                              upperSize=self.upper_size, 
                              upperProb=0.0, 
                              lowerProb=1.0)
    parents = selection(self.pop)
    self.valid_population(parents, 
                          self.selection_size, 
                          self.seq_length, 
                          survivorFlag=True)
    self.valid_couples(parents)
    self.sequences_from_bin(parents, self.lower_bin)
  


  def test_middle_bin_selection(self):
    selection = Overselection(numCouples=self.num_couples, 
                              upperSize=self.upper_size, 
                              upperProb=0.0, 
                              lowerProb=0.0)
    parents = selection(self.pop)
    self.valid_population(parents, 
                          self.selection_size, 
                          self.seq_length, 
                          survivorFlag=True)
    self.valid_couples(parents)
