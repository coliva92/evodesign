from abc import ABC
from Utils import EvoDesignTestCase
from evodesign.GA.Recombination.UniformCrossover import UniformCrossover
from evodesign.Population import Population





class RecombinationTests(EvoDesignTestCase, ABC):

  def setUp(self):
    self.sequences = [
      'AAAAAAAA',
      'CCCCCCCC',
      'DDDDDDDD',
      'EEEEEEEE',
      'FFFFFFFF',
      'GGGGGGGG',
      'HHHHHHHH',
      'IIIIIIII'
    ]
    self.parents = Population.create(self.sequences)
    self.num_parents = len(self.parents)
    self.seq_length = len(self.sequences[0])
    self.num_couples = self.num_parents // 2
  


  def valid_crossover(self, 
                      parents, 
                      children, 
                      numCouples):
    results = []
    for i in range(0, len(children), 2):
      mother = parents.iloc[i]['sequence']
      sister = children.iloc[i]['sequence']
      brother = children.iloc[i + 1]['sequence']
      x = [ a == b for a, b in zip(sister, mother) ]
      y = [ a == b for a, b in zip(brother, mother) ]
      z = [ a + b for a, b in zip(x, y) ]
      results.append(sum(z) == len(mother))
    num_correct_offspring_pairs = sum(results)
    self.assertEqual(num_correct_offspring_pairs, numCouples)





class UniformCrossoverTests(RecombinationTests):

  def test_uniform_crossover(self):
    recombination = UniformCrossover()
    children = recombination(self.parents)
    self.valid_population(children, self.num_parents, self.seq_length)
    self.valid_crossover(self.parents, children, self.num_couples)
