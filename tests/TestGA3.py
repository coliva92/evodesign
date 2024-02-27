from unittest import TestCase
from evodesign.Algorithms.GA3 import GA3
from evodesign.GA.Selection.Tournament import Tournament
from evodesign.GA.Recombination.UniformCrossover import UniformCrossover
from evodesign.GA.Mutation.RandomResetting import RandomResetting
from evodesign.Fitness.Rastrigin import Rastrigin
from evodesign.Prediction.Null import Null
import pandas as pd





class GA3Tests(TestCase):

  def setUp(self):
    self.algo = GA3(maxGenerations=100,
                    popSize=5,
                    predictor=Null(),
                    fitnessFn=Rastrigin(),
                    selection=Tournament(tournamentSize=2,
                                         sortColumns=[ 'fitness_rastrigin' ],
                                         sortAscending=[ False ]),
                    recombination=UniformCrossover(),
                    mutation=RandomResetting(),
                    sortColumns=[ 'fitness_rastrigin' ],
                    sortAscending=[ False ])
    self.pop = pd.DataFrame({
      'generation_id': 0,
      'sequence_id': [ 'A', 'B', 'C', 'D', 'E', 'X', 'Z' ],
      'fitness_rastrigin': [ -0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7 ],
      'survivor': [ True, True, True, True, True, False, False ]
    })
    self.children = pd.DataFrame({
      'generation_id': 1,
      'sequence_id': [ 'F', 'G', 'H', 'I', 'J' ],
      'fitness_rastrigin': [ -0.81, -0.82, -0.83, -0.84, -0.85 ],
      'survivor': [ False, False, False, False, False ]
    })
    self.correct_next_pop = pd.DataFrame({
      'generation_id': 1,
      'sequence_id': [ 'A', 'F', 'G', 'H', 'I', 'J' ],
      'fitness_rastrigin': [ -0.1, -0.81, -0.82, -0.83, -0.84, -0.85 ],
      'survivor': [ True, True, True, True, True, False ]
    })
    self.children_winner_first = pd.DataFrame({
      'generation_id': 1,
      'sequence_id': [ 'F', 'G', 'H', 'I', 'J' ],
      'fitness_rastrigin': [ -0.05, -0.82, -0.83, -0.84, -0.85 ],
      'survivor': [ False, False, False, False, False ]
    })
    self.correct_next_pop_winner_first = pd.DataFrame({
      'generation_id': 1,
      'sequence_id': [ 'F', 'G', 'H', 'I', 'J', 'A' ],
      'fitness_rastrigin': [ -0.05, -0.82, -0.83, -0.84, -0.85, -0.1 ],
      'survivor': [ True, True, True, True, True, False ]
    })
    self.children_winner_mid = pd.DataFrame({
      'generation_id': 1,
      'sequence_id': [ 'F', 'G', 'H', 'I', 'J' ],
      'fitness_rastrigin': [ -0.81, -0.82, -0.05, -0.84, -0.85 ],
      'survivor': [ False, False, False, False, False ]
    })
    self.correct_next_pop_winner_mid = pd.DataFrame({
      'generation_id': 1,
      'sequence_id': [ 'H', 'F', 'G', 'I', 'J', 'A' ],
      'fitness_rastrigin': [ -0.05, -0.81, -0.82, -0.84, -0.85, -0.1 ],
      'survivor': [ True, True, True, True, True, False ]
    })
    self.children_winner_last = pd.DataFrame({
      'generation_id': 1,
      'sequence_id': [ 'F', 'G', 'H', 'I', 'J' ],
      'fitness_rastrigin': [ -0.81, -0.82, -0.83, -0.84, -0.05 ],
      'survivor': [ False, False, False, False, False ]
    })
    self.correct_next_pop_winner_last = pd.DataFrame({
      'generation_id': 1,
      'sequence_id': [ 'J', 'F', 'G', 'H', 'I', 'A' ],
      'fitness_rastrigin': [ -0.05, -0.81, -0.82, -0.83, -0.84, -0.1 ],
      'survivor': [ True, True, True, True, True, False ]
    })
  


  def test_replacement(self):
    next_pop = self.algo.replacement(self.pop, self.children, [])
    comparison = next_pop == self.correct_next_pop
    rows_equal = comparison.all(axis=1)
    self.assertEqual(rows_equal.all(), True)
  


  def test_replacement_winner_first(self):
    next_pop = self.algo.replacement(self.pop, self.children_winner_first, [])
    comparison = next_pop == self.correct_next_pop_winner_first
    rows_equal = comparison.all(axis=1)
    self.assertEqual(rows_equal.all(), True)



  def test_replacement_winner_mid(self):
    next_pop = self.algo.replacement(self.pop, self.children_winner_mid, [])
    comparison = next_pop == self.correct_next_pop_winner_mid
    rows_equal = comparison.all(axis=1)
    self.assertEqual(rows_equal.all(), True)
  


  def test_replacement_winner_last(self):
    next_pop = self.algo.replacement(self.pop, self.children_winner_last, [])
    comparison = next_pop == self.correct_next_pop_winner_last
    rows_equal = comparison.all(axis=1)
    self.assertEqual(rows_equal.all(), True)
