from unittest import TestCase
from evodesign.Chain import Chain





class ChainTests(TestCase):

  def setUp(self):
    self.structure = Chain.load_structure('tests/8hjc.pdb')
    correct_sequence = 'CLESGTSCIPGAQHNCCSGVCVPIVTIFYGVCY'
    self.correct_length = len(correct_sequence)
    self.correct_backbone = self.correct_length * list(Chain.BACKBONE_ATOMS)



  def test_backbone_atoms_retrieval(self):
    backbone = Chain.backbone_atoms(self.structure)
    retrieved_backbone = [ atom.get_name() for atom in backbone ]
    self.assertEqual(retrieved_backbone, self.correct_backbone)



  def test_backbone_length_retrieval(self):
    seq_length = Chain.length(self.structure)
    self.assertIsInstance(seq_length, int)
    self.assertEqual(seq_length, self.correct_length)
