from evodesign.GA.Recombination.UniformCrossover import UniformCrossover
from evodesign.Population import Population

def valid_crossover(parents, children, numCouples):
  results = []
  for i in range(0, len(children), 2):
    mother = parents.iloc[i]['sequence']
    sister = children.iloc[i]['sequence']
    brother = children.iloc[i + 1]['sequence']
    x = [ a == b for a, b in zip(sister, mother) ]
    y = [ a == b for a, b in zip(brother, mother) ]
    z = [ a + b for a, b in zip(x, y) ]
    results.append(sum(z) == len(mother))
  return sum(results) == numCouples

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
parents = Population.create(sequences)
recombination = UniformCrossover()
children = recombination(parents)
test_passed = valid_crossover(parents, children, len(sequences) // 2)
print('PASSED' if test_passed else 'FAILED')
