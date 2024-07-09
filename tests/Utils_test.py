from unittest import TestCase
import evodesign.Utils as op
import pandas as pd
import numpy as np





class UtilsTest(TestCase):

  def setUp(self):
    self.array_A = np.array([ 1.0, 1.5, 2.0, 0.0 ])
    self.array_B = np.array([ 1.0, 1.5, 1.5, 0.0 ])
    self.array_C = np.array([ 1.0, 1.5, 1.5, 2.0 ])
    self.series_A = pd.Series({ 'value1': 1.0, 'value2': 1.0 })
    self.series_B = pd.Series({ 'value1': 0.5, 'value2': 1.0 })
    self.series_C = pd.Series({ 'value1': 1.0, 'value2': 0.5 })
    self.pop_A = pd.DataFrame({
      'value1': [ 1.0, 1.0, 0.8, 0.7, 0.6 ],
      'value2': [ 0.1, 0.5, 0.1, 0.4, 0.2 ]
    })
    self.pop_B = pd.DataFrame({
      'value1': [ 1.0, 0.8, 0.7, 0.7, 0.7 ],
      'value2': [ 0.7, 0.1, 0.2, 0.3, 0.5 ]
    })
    self.sorted_pop = pd.DataFrame({
      'value1': [ 1.0, 1.0, 1.0, 0.8, 0.8, 0.7, 0.7, 0.7, 0.7, 0.6 ],
      'value2': [ 0.1, 0.5, 0.7, 0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.2 ]
    })
    self.sorted_pop2 = pd.DataFrame({
      'value1': [ 1.0, 1.0, 1.0, 0.8, 0.8, 0.7, 0.7, 0.7 ],
      'value2': [ 0.1, 0.5, 0.7, 0.1, 0.1, 0.2, 0.3, 0.5 ]
    })
    self.sorted_pop3 = pd.DataFrame({
      'value1': [ 1.0, 1.0, 1.0, 0.8, 0.8, 0.7, 0.7, 0.6 ],
      'value2': [ 0.1, 0.5, 0.7, 0.1, 0.1, 0.2, 0.4, 0.2 ]
    })



  def test_one_dominates_other(self):
    A_dominates_B = op.dominates(self.array_A, self.array_B)
    self.assertEqual(A_dominates_B, True)
  


  def test_other_not_dominates_one(self):
    B_dominates_A = op.dominates(self.array_B, self.array_A)
    self.assertEqual(B_dominates_A, False)
  


  def test_one_not_dominates_one(self):
    A_dominates_A = op.dominates(self.array_A, self.array_A)
    self.assertEqual(A_dominates_A, False)
  


  def test_one_not_dominates_third(self):
    A_dominates_C = op.dominates(self.array_A, self.array_C)
    self.assertEqual(A_dominates_C, False)
  

  def test_third_not_dominates_one(self):
    C_dominates_A = op.dominates(self.array_C, self.array_A)
    self.assertEqual(C_dominates_A, False)
  


  def test_first_before_second_value1(self):
    A_before_B = op.is_sorted_before(self.series_A, 
                                self.series_B,
                                [ 'value1' ],
                                [ False ])
    self.assertEqual(A_before_B, True)
  


  def test_second_not_before_first_value1(self):
    B_before_A = op.is_sorted_before(self.series_B,
                                self.series_A,
                                [ 'value1' ],
                                [ False ])
    self.assertEqual(B_before_A, False)
  


  def test_first_before_third_value1(self):
    A_before_C = op.is_sorted_before(self.series_A,
                                self.series_C,
                                [ 'value1' ],
                                [ False ])
    self.assertEqual(A_before_C, True)
  


  def test_first_not_before_third_value1_value2(self):
    A_before_C = op.is_sorted_before(self.series_A, 
                                self.series_C,
                                [ 'value1', 'value2' ],
                                [ False, True ])
    self.assertEqual(A_before_C, False)



  def test_third_before_first_value1_value2(self):
    C_before_A = op.is_sorted_before(self.series_C, 
                                self.series_A,
                                [ 'value1', 'value2' ],
                                [ False, True ])
    self.assertEqual(C_before_A, True)



  def test_first_before_first_value1_value2(self):
    A_before_A = op.is_sorted_before(self.series_A,
                                self.series_A,
                                [ 'value1', 'value2' ],
                                [ False, True ])
    self.assertEqual(A_before_A, True)
  


  def test_merge(self):
    merged = op.merge(self.pop_A, 
                      self.pop_B,
                      [ 'value1', 'value2' ],
                      [ False, True ])
    df_equal = self.sorted_pop.equals(merged)
    self.assertEqual(df_equal, True)
  


  def test_merge_fixed_size(self):
    size = 5
    merged = op.merge(self.pop_A,
                      self.pop_B,
                      [ 'value1', 'value2' ],
                      [ False, True ],
                      size)
    df_equal = self.sorted_pop.iloc[:size].equals(merged)
    self.assertEqual(df_equal, True)
  


  def test_merged_first_smaller(self):
    merged = op.merge(self.pop_A.iloc[:3],
                      self.pop_B,
                      [ 'value1', 'value2' ],
                      [ False, True ])
    df_equal = self.sorted_pop2.equals(merged)
    self.assertEqual(df_equal, True)
  


  def test_merged_second_smaller(self):
    merged = op.merge(self.pop_A,
                      self.pop_B.iloc[:3],
                      [ 'value1', 'value2' ],
                      [ False, True ])
    df_equal = self.sorted_pop3.equals(merged)
    self.assertEqual(df_equal, True)
  


  def test_merged_first_smaller_fixed_size(self):
    size = 5
    merged = op.merge(self.pop_A.iloc[:3],
                      self.pop_B,
                      [ 'value1', 'value2' ],
                      [ False, True ],
                      size)
    df_equal = self.sorted_pop2.iloc[:size].equals(merged)
    self.assertEqual(df_equal, True)
  


  def test_merged_second_smaller_fixed_size(self):
    size = 5
    merged = op.merge(self.pop_A,
                      self.pop_B.iloc[:3],
                      [ 'value1', 'value2' ],
                      [ False, True ],
                      size)
    df_equal = self.sorted_pop3.iloc[:size].equals(merged)
    self.assertEqual(df_equal, True)
