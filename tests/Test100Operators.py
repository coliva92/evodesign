from unittest import TestCase
import evodesign.Operators as op
import numpy as np





class OperatorsTests(TestCase):

  def setUp(self):
    self.array_A = np.array([ 1.0, 1.5, 2.0, 0.0 ])
    self.array_B = np.array([ 1.0, 1.5, 1.5, 0.0 ])
    self.array_C = np.array([ 1.0, 1.5, 1.5, 2.0 ])



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
