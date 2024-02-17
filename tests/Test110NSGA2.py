from Utils import EvoDesignTestCase
from evodesign.Algorithms.NSGA2 import NSGA2
from evodesign.Prediction.Null import Null
from evodesign.GA.Recombination.UniformCrossover import UniformCrossover
from evodesign.GA.Mutation.RandomResetting import RandomResetting
from evodesign.Fitness.Gdt import Gdt
from evodesign.Fitness.Rmsd import Rmsd
from evodesign.Fitness.Rastrigin import Rastrigin
import pandas as pd





class NSGA2Tests(EvoDesignTestCase):

  def setUp(self):
    self.algo = NSGA2(
      maxGenerations=100,
      popSize=10,
      predictor=Null(),
      fitnessFns=[ Rmsd(), Gdt(), Rastrigin() ],
      recombination=UniformCrossover(),
      mutation=RandomResetting()
    )



  def test_sort_by_rank_and_distance(self):
    unsorted_data = {
      'rank': [ 2, 1, 3, 2, 1 ],
      'distance': [ 1.0, 1.0, 3.0, 2.0, 2.0 ]
    }
    sorted_data = {
      'rank': [ 1, 1, 2, 2, 3 ],
      'distance': [ 2.0, 1.0, 2.0, 1.0, 3.0 ]
    }
    unsorted_df = pd.DataFrame(unsorted_data)
    sorted_df = pd.DataFrame(sorted_data)
    df = unsorted_df.sort_values(by=[ 'rank', 'distance' ],
                                 ascending=[ True, False ],
                                 ignore_index=True)
    comparison = df == sorted_df
    rows_equal = comparison.all(axis=1)
    self.assertTrue(rows_equal.all())
  


  def test_non_domination_rank(self):
    unranked_data = {
      'generation_id': 9 * [ 0 ],
      'sequence_id': [ f'prot_0000_000{i}' for i in range(9) ],
      'sequence': [ 
        'AAAAAAAA', 'CCCCCCCC', 'DDDDDDDD', 'EEEEEEEE', 'FFFFFFFF', 'GGGGGGGG', 
        'HHHHHHHH', 'IIIIIIII', 'KKKKKKKK'
      ],
      'survivor': 9 * [ True ],
      'fitness_rmsd': [
        -1.6, -1.8, -1.5, -1.5, -2.5, -1.5, -1.5, -1.8, -1.5
      ],
      'fitness_gdt': [
        0.71, 0.71, 0.8, 0.71, 0.67, 0.79, 0.74, 0.71, 0.78
      ],
      'fitness_rastrigin': [
        -30, -36, -28, -31, -33, -25, -33, -29, -24
      ]
    }
    ranked_data = {
      'generation_id': 9 * [ 0 ],
      'sequence_id': [ f'prot_0000_000{i}' for i in range(9) ],
      'sequence': [ 
        'AAAAAAAA', 'CCCCCCCC', 'DDDDDDDD', 'EEEEEEEE', 'FFFFFFFF', 'GGGGGGGG', 
        'HHHHHHHH', 'IIIIIIII', 'KKKKKKKK',
      ],
      'survivor': 9 * [ True ],
      'fitness_rmsd': [
        -1.6, -1.8, -1.5, -1.5, -2.5, -1.5, -1.5, -1.8, -1.5
      ],
      'fitness_gdt': [
        0.71, 0.71, 0.8, 0.71, 0.67, 0.79, 0.74, 0.71, 0.78
      ],
      'fitness_rastrigin': [
        -30, -36, -28, -31, -33, -25, -33, -29, -24
      ],
      'rank': [ 2, 3, 1, 2, 3, 1, 2, 2, 1 ]
    }
    unranked_df = pd.DataFrame(unranked_data)
    ranked_df = pd.DataFrame(ranked_data)
    correct_fronts = [
      ranked_df.loc[[ 2, 5, 8 ]], 
      ranked_df.loc[[ 0, 3, 6, 7 ]],
      ranked_df.loc[[ 1, 4 ]]
    ]
    pop, fronts = self.algo.non_domination_rank(unranked_df)
    self.assertEqual(len(fronts), 3)
    self.assertIn('rank', pop.columns)
    comparison = pop['rank'] == ranked_df['rank']
    self.assertTrue(comparison.all())
    for i in range(3):
      front, correct = fronts[i], correct_fronts[i]
      comparison = front == correct
      rows_equal = comparison.all(axis=1)
      self.assertTrue(rows_equal.all())
