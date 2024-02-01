from evodesign.GA.Selection.Overselection import Overselection
from evodesign.Population import Population
import pandas as pd

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
test1_passed = valid_couples(parents, num_couples, upper_bin)

selection = Overselection(numSelectedCouples=num_couples, 
                          upperSize=upper_size, 
                          upperProb=0., 
                          lowerProb=1.)
parents = selection(pop)
test2_passed = valid_couples(parents, num_couples, lower_bin)

selection = Overselection(numSelectedCouples=num_couples, 
                          upperSize=upper_size, 
                          upperProb=0., 
                          lowerProb=0.)
parents = selection(pop)
test3_passed = valid_couples(parents, num_couples, sequences)

all_tests_passed = test1_passed and \
                   test2_passed and \
                   test3_passed
print('PASSED' if all_tests_passed else 'FAILED')
