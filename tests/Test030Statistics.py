from unittest import TestCase
from evodesign.Statistics import Statistics
from evodesign.Population import Population





class StatisticsTests(TestCase):

  def setUp(self):
    self.pop = Population.create([
      'AAAAA',
      'CAAAC',
      'DAACC',
      'EACCE',
      'FCCDE'
    ])
  


  def test_amino_acid_loss(self):
    x = Statistics.average_amino_acid_loss(self.pop)
    self.assertEqual(x, 17.0)
  


  def test_sequence_identity(self):
    # 3 2 1 0 => 1.5 * (4/(4+3+2+1)) => 0.6
    #   3 1 0 => 1.3333 * (3/(4+3+2+1)) => 0.4
    #     2 0 => 1.0 * (2/(4+3+2+1)) => 0.2
    #       2 => 2.0 * (1/(4+3+2+1)) => 0.2
    # 1.4
    x = Statistics.average_sequence_identity(self.pop)
    self.assertEqual(x, 1.4)
