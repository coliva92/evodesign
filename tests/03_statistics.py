from evodesign.Statistics import Statistics
from evodesign.Population import Population
import pandas as pd

sequences = [
  'AAAAA',
  'CAAAC',
  'DAACC',
  'EACCE',
  'FCCDE'
]
pop = Population.create(sequences)
a = Statistics.average_missing_amino_acids(pop)
test1_passed = a == 17.
b = Statistics.average_sequence_identity(pop)
test2_passed = b == 1.4
print('PASSED' if test1_passed and test2_passed else 'FAILED')

# 3 2 1 0 => 1.5 * (4/(4+3+2+1)) => 0.6
#   3 1 0 => 1.3333 * (3/(4+3+2+1)) => 0.4
#     2 0 => 1.0 * (2/(4+3+2+1)) => 0.2
#       2 => 2.0 * (1/(4+3+2+1)) => 0.2
# 1.4
