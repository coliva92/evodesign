from Utils import EvoDesignTestCase
import evodesign.Population as Population





class PopulationTests(EvoDesignTestCase):

  def setUp(self):
    self.correct_sequences = [
      'AAAAAA',
      'CCCCCC',
      'DDDDDD',
      'EEEEEE'
    ]
  


  def test_pop_creation(self):
    pop = Population.create(self.correct_sequences)
    self.valid_population(pop, 
                          len(self.correct_sequences),
                          len(self.correct_sequences[0]))
    pop_sequences = [ s for _, s in pop['sequence'].items() ]
    self.assertEqual(pop_sequences, self.correct_sequences)
  


  def test_random_pop_creation(self):
    pop_size, seq_length = 5, 10
    pop = Population.create_random(pop_size, seq_length)
    self.valid_population(pop, pop_size, seq_length)
