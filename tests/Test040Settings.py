from unittest import TestCase
from evodesign.Settings import Settings
from evodesign.Fitness.Rmsd import Rmsd
from evodesign.Fitness.Gdt import Gdt
from evodesign.Fitness.Experimental.Rastrigin import Rastrigin
from evodesign.Fitness.Experimental.RmsdCyclization import RmsdCyclization
from evodesign.Fitness.Experimental.SideChainPacking import SideChainPacking
from evodesign.GA.Selection.Tournament import Tournament
from evodesign.Workspace import Workspace
import math






class SettingsTests(TestCase):

  def setUp(self):
    Workspace('foo', 'bar')
  


  def test_fitness_rmsd_settings_retrieval(self):
    correct_settings = {
      'Fitness.Rmsd': {
        'upperBound': -0.5
      }
    }
    rmsd = Rmsd()
    settings = rmsd.settings()
    self.assertEqual(settings, correct_settings)
  


  def test_fitness_gdt_settings_retrieval(self):
    correct_settings = {
      'Fitness.Gdt': {
        'upperBound': 0.95,
        'cutoffs': [ 1.0, 2.0, 4.0, 8.0 ]
      }
    }
    gdt = Gdt()
    settings = gdt.settings()
    self.assertEqual(settings, correct_settings)
  


  def test_fitness_rastrigin_settings_retrieval(self):
    correct_settings = {
      'Fitness.Rastrigin': {
        'upperBound': 0.0,
        'windowWidth': 3
      }
    }
    rastrigin = Rastrigin()
    settings = rastrigin.settings()
    self.assertEqual(settings, correct_settings)
  


  def test_fitness_cyclization_settings_retrieval(self):
    correct_settings = {
      'Fitness.Experimental.RmsdCyclization': {
        'upperBound': 1.0,
        'rmsdWeight': 1.0,
        'cycWeight': 1.0
      }
    }
    cyclization = RmsdCyclization()
    settings = cyclization.settings()
    self.assertEqual(settings, correct_settings)
  


  def test_fitness_sc_packing_settings_retrieval(self):
    correct_settings = {
      'Fitness.Experimental.SideChainPacking': {
        'upperBound': -math.inf,
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
          'Fitness.Gdt': {
            'cutoffs': [ 0.5, 1.0 ]
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
          'GA.Recombination.MiddlePointCrossover': {
            'probability': 1.0
          }
        },
        'mutation': {
          'GA.Mutation.Switch': {}
        }
      }
    }
    algo = Settings.parse(correct_settings)
    correct_settings['Algorithms.GA2']['betterOffspringBias'] = 1.0
    correct_settings['Algorithms.GA2']['sortColumns'] = [ 'fitness_gdt', 'rmsd', 'plddt' ]
    correct_settings['Algorithms.GA2']['sortAscending'] = [ False, True, False ]
    correct_settings['Algorithms.GA2']['fitnessFn']['Fitness.Gdt']['upperBound'] = 0.95
    correct_settings['Algorithms.GA2']['predictor']['Prediction.AlphaFold']['maxTemplateDate'] = '2020-05-14'
    correct_settings['Algorithms.GA2']['predictor']['Prediction.AlphaFold']['modelPreset'] = 'monomer'
    correct_settings['Algorithms.GA2']['predictor']['Prediction.AlphaFold']['dbPreset'] = 'reduced_dbs'
    correct_settings['Algorithms.GA2']['mutation']['GA.Mutation.Switch']['mutProb'] = 1.0
    correct_settings['Algorithms.GA2']['mutation']['GA.Mutation.Switch']['numSwitches'] = 1
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
