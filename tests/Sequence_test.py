from Utils import EvoDesignTestCase
import evodesign.Sequence as Sequence
import evodesign.Random as Random





class SequenceTests(EvoDesignTestCase):

  def test_random_sequence_generation(self):
    sequence = Sequence.create_random(10)
    self.valid_sequence(sequence, 10)



  def test_amino_acid_exchange(self):
    rng = Random.generator()
    old_letter = rng.choice(list(Sequence.AMINO_ACIDS))
    new_letter = Sequence.swap_letter(old_letter)
    self.assertIsInstance(new_letter, str)
    new_letter_length = len(new_letter)
    self.assertEqual(new_letter_length, 1)
    self.assertNotEqual(new_letter, old_letter)
    new_letter_is_amino_acid = new_letter in Sequence.AMINO_ACIDS
    self.assertEqual(new_letter_is_amino_acid, True)
