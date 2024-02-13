from abc import ABC
from unittest import TestCase
from evodesign.Population import Population
from evodesign.Sequence import Sequence
import pandas as pd





class EvoDesignTestCase(TestCase, ABC):

  def valid_generation_id(self, 
                          generationId, 
                          value):
    self.assertEqual(generationId, value)
    return True



  def valid_sequence_id(self, sequenceId):
    self.assertIsInstance(sequenceId, str)
    self.assertRegex(sequenceId, r'^prot_\d{4}_\d{4}$')
    return True



  def valid_sequence(self,
                     sequence, 
                     length):
    self.assertIsInstance(sequence, str)
    sequence_length = len(sequence)
    self.assertEqual(sequence_length, length)
    num_amino_acids = sum([ i in Sequence.AMINO_ACIDS for i in sequence ])
    self.assertEqual(num_amino_acids, length)
    return True



  def valid_survivor_flag(self, 
                          survivor, 
                          value = False):
    self.assertEqual(survivor, value)
    return True



  def valid_population(self, 
                       population, 
                       popSize, 
                       seqLength,
                       generationId = 0,
                       survivorFlag = False):
    self.assertIsInstance(population, pd.DataFrame)
    size = len(population)
    self.assertEqual(size, popSize)
    pop_has_sequence = 'sequence' in population
    pop_has_generation_id = 'generation_id' in population
    pop_has_survivor = 'survivor' in population
    pop_has_sequence_id = 'sequence_id' in population
    self.assertEqual(pop_has_sequence, True)
    self.assertEqual(pop_has_generation_id, True)
    self.assertEqual(pop_has_survivor, True)
    self.assertEqual(pop_has_sequence_id, True)
    num_valid_sequences = \
      population['sequence'].apply(self.valid_sequence, length=seqLength).sum()
    num_valid_generation_ids = \
      population['generation_id'].apply(self.valid_generation_id, 
                                        value=generationId).sum()
    num_valid_survivor_flags = \
      population['survivor'].apply(self.valid_survivor_flag, 
                                   value=survivorFlag).sum()
    num_valid_sequence_ids = sum([
      self.valid_sequence_id(s)
      for _, s in population['sequence_id'].items()
    ])
    self.assertEqual(num_valid_sequences, popSize)
    self.assertEqual(num_valid_generation_ids, popSize)
    self.assertEqual(num_valid_survivor_flags, popSize)
    self.assertEqual(num_valid_sequence_ids, popSize)
    return True
