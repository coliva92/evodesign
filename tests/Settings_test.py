from unittest import TestCase
import evodesign.Settings as Settings
from evodesign.Fitness.LinearCombination import LinearCombination
from evodesign.Metrics.Gdt import Gdt
from evodesign.Metrics.Rmsd import Rmsd
from evodesign.Fitness.Experimental.Rastrigin import Rastrigin
from evodesign.Metrics.Experimental.SideChainPacking import SideChainPacking
from evodesign.GA.Selection.Tournament import Tournament
from evodesign.Workspace import Workspace
import math






class SettingsTests(TestCase):

  def setUp(self):
    Workspace('foo', 'bar')
  


  def test_fitness_linear_combination_settings_retrieval(self):
    correct_settings = {
      'Fitness.LinearCombination': {
        'upperBound': 0.99,
        'coefficients': [ 1.0 ],
        'terms': [
          {
            "Metrics.Gdt": {
              "cutoffs": [ 1.0, 2.0, 4.0, 8.0 ],
              "rmsdCalculator": {
                "Metrics.Rmsd": {}
              }
            }
          }
        ]
      }
    }
    fitness_fn = LinearCombination(upperBound=0.99, terms=[ Gdt(rmsdCalculator=Rmsd()) ])
    settings = fitness_fn.settings()
    self.assertEqual(settings, correct_settings)
  


  def test_rastrigin_settings_retrieval(self):
    correct_settings = {
      'Fitness.Experimental.Rastrigin': {
        'upperBound': 0.0,
        'windowWidth': 3,
        'terms': []
      }
    }
    rastrigin = Rastrigin()
    settings = rastrigin.settings()
    self.assertEqual(settings, correct_settings)
  


  def test_sc_packing_settings_retrieval(self):
    correct_settings = {
      'Metrics.Experimental.SideChainPacking': {
        'scwrlExecutablePath': './scwrl4/Scwrl4'
      }
    }
    sc_packing = SideChainPacking()
    settings = sc_packing.settings()
    self.assertEqual(settings, correct_settings)
  


  def test_parse_algorithm(self):
    correct_settings = {
      'Algorithms.GA2': {
        'maxGenerations': 1000,
        'popSize': 100,
        'elitismSize': 20,
        'sortColumns': [ 'fitness_gdt', 'rmsd', 'plddt' ],
        'sortAscending': [ False, True, False ],
        'predictor': {
          'Prediction.AlphaFold': {
            'fakeMsaScript': 'fakemsa.py',
            'alphafoldScript': 'run_docker.py',
            'mgnifyDbPath': '/media/my_usb/mgnify',
            'dataDir': '/home/com/my_usb'
          }
        },
        'fitnessFn': {
          "Fitness.LinearCombination": {
            "upperBound": 0.95,
            "terms": [
              {
                'Metrics.Gdt': {
                  'cutoffs': [ 0.5, 1.0 ],
                  'rmsdCalculator': { 'Metrics.Rmsd': {} }
                }
              }
            ]
          }
          
        },
        'selection': {
          'GA.Selection.Tournament': {
            'tournamentSize': 5,
            'sortColumns': [ 'rank', 'distance' ],
            'sortAscending': [ True, False ]
          }
        },
        'recombination': {
          'GA.Recombination.MiddlePointCrossover': {}
        },
        'mutation': {
          'GA.Mutation.Swap': {}
        }
      }
    }
    algo = Settings.parse(correct_settings)
    correct_settings['Algorithms.GA2']['betterOffspringBias'] = 1.0
    correct_settings['Algorithms.GA2']['sortColumns'] = [ 'fitness_gdt', 'rmsd', 'plddt' ]
    correct_settings['Algorithms.GA2']['sortAscending'] = [ False, True, False ]
    correct_settings['Algorithms.GA2']['fitnessFn']['Fitness.LinearCombination']['upperBound'] = 0.95
    correct_settings['Algorithms.GA2']['predictor']['Prediction.AlphaFold']['maxTemplateDate'] = '2020-05-14'
    correct_settings['Algorithms.GA2']['predictor']['Prediction.AlphaFold']['modelPreset'] = 'monomer'
    correct_settings['Algorithms.GA2']['predictor']['Prediction.AlphaFold']['dbPreset'] = 'reduced_dbs'
    correct_settings['Algorithms.GA2']['mutation']['GA.Mutation.Swap']['mutProb'] = 1.0
    correct_settings['Algorithms.GA2']['mutation']['GA.Mutation.Swap']['numSwaps'] = 1
    correct_settings['Algorithms.GA2']['selection']['GA.Selection.Tournament']['elitism'] = False
    settings = algo.settings()
    self.assertEqual(settings, correct_settings)
  


  def test_tournament_selection_settings(self):
    correct_settings = {
      'GA.Selection.Tournament': {
        'tournamentSize': 3,
        'sortColumns': [ 'rank', 'distance' ],
        'sortAscending': [ True, False ],
        'elitism': True
      }
    }
    tournament = Tournament(elitism=True,
                            tournamentSize=3,
                            sortColumns=[ 'rank', 'distance' ],
                            sortAscending=[ True, False ])
    settings = tournament.settings()
    self.assertEqual(settings, correct_settings)
