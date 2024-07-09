from unittest import TestCase
from evodesign.Prediction.ESMFoldRemoteApi import ESMFoldRemoteApi
import evodesign.Chain as Chain
import os





class PredictionTests(TestCase):

  def setUp(self):
    self.sequence = 'GREETINGSTRAVELER'
    self.filepath = 'tests/esmfold_prediction.pdb'
    self.predictor = ESMFoldRemoteApi()
    self.correct_backbone = len(self.sequence) * list(Chain.BACKBONE_ATOMS)
  


  def tearDown(self):
    os.remove(self.filepath)
  


  def test_prediction(self):
    backbone, plddt = self.predictor(self.sequence, self.filepath)
    self.assertGreaterEqual(plddt, 0.0)
    self.assertLessEqual(plddt, 1.0)
    backbone_atoms = [ atom.get_name() for atom in backbone ]
    self.assertEqual(backbone_atoms, self.correct_backbone)
    self.assertTrue(os.path.isfile(self.filepath))
