from abc import ABC
from Utils import EvoDesignTestCase
from evodesign.GA.Mutation.RandomResetting import RandomResetting
from evodesign.Population import Population





class MutationTests(EvoDesignTestCase, ABC):

  def setUp(self):
    self.sequences = [
      'AAAAAAAA',
      'CCCCCCCC',
      'DDDDDDDD'
    ]
    self.original = Population.create(self.sequences)
    self.pop_size = len(self.original)
    self.seq_length = len(self.sequences[0])





class RandomResettingTests(MutationTests):

  def test_random_resetting(self):
    mutation = RandomResetting(exchangeProb=0.5)
    mutated = mutation(self.original.copy())
    self.valid_population(mutated, self.pop_size, self.seq_length)
    mutated_size = len(mutated)
    self.assertEqual(mutated_size, self.pop_size)
    num_mutated_sequences = (mutated['sequence'] != self.original['sequence']).sum()
    self.assertEqual(num_mutated_sequences, self.pop_size)
    num_equal_sequence_ids = (mutated['sequence_id'] == self.original['sequence_id']).sum()
    self.assertEqual(num_equal_sequence_ids, self.pop_size)
