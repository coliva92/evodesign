from evodesign.GA.Selection.Overselection import Overselection
from evodesign.Population import Population
from evodesign.Sequence import Sequence
import pandas as pd
import numpy as np

def valid_generation_id(generationId, value):
  return (type(generationId) == int or type(generationId) == np.int64) and \
         generationId == value

def valid_sequence_id(sequenceId, generationId):
  return type(sequenceId) == str and \
         len(sequenceId.split('_')) == 3 and \
         'prot' in sequenceId and \
         len(sequenceId) == 14 and \
         int(float(sequenceId[5:9])) == generationId

def valid_sequence(sequence, length):
  return type(sequence) == str and \
         len(sequence) == length and \
         sum([ i in Sequence.AMINO_ACIDS for i in sequence ]) == len(sequence)

def valid_survivor_flag(survivor):
  return (type(survivor) == bool or type(survivor) == np.bool_) and \
         survivor == True

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
           valid_sequence_id(s, 0) 
           for _, s in population['sequence_id'].items() 
         ]) == n

def couples_of_different_sequences(parents, numCouples):
  return sum([ 
    parents.iloc[i]['sequence'] != parents.iloc[i + 1]['sequence'] 
    for i in range(0, len(parents), 2) 
  ]) == numCouples

def valid_couples(parents, numCouples, bin):
  n = 2 * numCouples
  return isinstance(parents, pd.DataFrame) and \
         len(parents) == n and \
         sum([ row['sequence'] in bin for _, row in parents.iterrows() ]) == n and \
         couples_of_different_sequences(parents, numCouples)

sequences = [
  'AAAAAAAA',
  'CCCCCCCC',
  'DDDDDDDD',
  'EEEEEEEE',
  'FFFFFFFF',
  'GGGGGGGG',
  'HHHHHHHH',
  'IIIIIIII'
]
num_couples = 3
upper_size = 3
upper_bin = { seq for seq in sequences[:upper_size] }
lower_bin = { seq for seq in sequences[upper_size:] }
pop = Population.create(sequences)
pop['survivor'] = True
selection = Overselection(numSelectedCouples=num_couples, 
                          upperSize=upper_size, 
                          upperProb=1., 
                          lowerProb=0.)
parents = selection(pop)
test1_passed = valid_population(parents, 2 * num_couples, len(sequences[0])) and \
               valid_couples(parents, num_couples, upper_bin)

selection = Overselection(numSelectedCouples=num_couples, 
                          upperSize=upper_size, 
                          upperProb=0., 
                          lowerProb=1.)
parents = selection(pop)
test2_passed = valid_population(parents, 2 * num_couples, len(sequences[0])) and \
               valid_couples(parents, num_couples, lower_bin)

selection = Overselection(numSelectedCouples=num_couples, 
                          upperSize=upper_size, 
                          upperProb=0., 
                          lowerProb=0.)
parents = selection(pop)
test3_passed = valid_population(parents, 2 * num_couples, len(sequences[0])) and \
               valid_couples(parents, num_couples, sequences)

all_tests_passed = test1_passed and \
                   test2_passed and \
                   test3_passed
print('PASSED' if all_tests_passed else 'FAILED')
