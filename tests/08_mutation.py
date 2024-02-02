from evodesign.GA.Mutation.RandomResetting import RandomResetting
from evodesign.Population import Population
from evodesign.Sequence import Sequence
import pandas as pd

def valid_generation_id(generationId, value):
  return type(generationId) == int and \
         generationId == value

def valid_sequence_id(sequenceId, generationId, rowIdx):
  return type(sequenceId) == str and \
         len(sequenceId.split('_')) == 3 and \
         'prot' in sequenceId and \
         len(sequenceId) == 14 and \
         int(float(sequenceId[5:9])) == generationId and \
         int(float(sequenceId[10:])) == rowIdx

def valid_sequence(sequence, length):
  return type(sequence) == str and \
         len(sequence) == length and \
         sum([ i in Sequence.AMINO_ACIDS for i in sequence ]) == len(sequence)

def valid_survivor_flag(survivor):
  return type(survivor) == bool and \
         survivor == False

def valid_population(population, popSize, seqLen):
  n = len(population)
  return isinstance(population, pd.DataFrame) and \
         len(population) == popSize and \
         'sequence' in population and \
         'generation_id' in population and \
         'survivor' in population and \
         'sequence_id' in population and \
         population['sequence'].apply(valid_sequence, length=seqLen).sum() == n and \
         population['generation_id'].apply(valid_generation_id, value=0).sum() == n and \
         population['survivor'].apply(valid_survivor_flag).sum() == n and \
         sum([ 
           valid_sequence_id(s, 0, i) 
           for i, s in population['sequence_id'].items() 
         ]) == n

sequences = [
  'AAAAAAAA',
  'CCCCCCCC',
  'DDDDDDDD'
]
before = Population.create(sequences)
mutation = RandomResetting(exchangeProb=0.5)
after = mutation(before.copy())
test_passed = valid_population(after, len(sequences), len(sequences[0])) and \
              len(after) == len(before) and \
              (after['sequence_id'] == before['sequence_id']).sum() == len(before) and \
              (after['sequence'] != before['sequence']).sum() == len(before)
print('PASSED' if test_passed else 'FAILED')
