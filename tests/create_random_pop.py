import evodesign.Population as Population
import pandas as pd
import evodesign.Sequence as Sequence


def valid_sequence(sequence):
  return type(sequence) == str and \
         len(sequence) == 10 and \
         sum([ i in Sequence.AMINO_ACIDS for i in sequence ]) == len(sequence)

pop = Population.create_random(5, 10)
test_passed = isinstance(pop, pd.DataFrame) and \
              len(pop.index) == 5 and \
              'Sequence' in pop and \
              pop['Sequence'].apply(valid_sequence).sum() == len(pop.index)
print('PASSED' if test_passed else 'FAILED')
