import pandas as pd
from evodesign.Population import Population
from evodesign.Sequence import Sequence


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
  'AAAAAA',
  'CCCCCC',
  'DDDDDD',
  'EEEEEE'
]
pop = Population.create(sequences)
pop_seqs = { s for _, s in pop['sequence'].items() }
test1_passed = valid_population(pop, 4, 6) and \
               'AAAAAA' in pop_seqs and \
               'CCCCCC' in pop_seqs and \
               'DDDDDD' in pop_seqs and \
               'EEEEEE' in pop_seqs
pop = Population.create_random(5, 10)
test2_passed = valid_population(pop, 5, 10)
print('PASSED' if test1_passed and test2_passed else 'FAILED')
