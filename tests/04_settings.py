from evodesign.Settings import Settings
from evodesign.Algorithms.PGPD import PGPD
from evodesign.Prediction.AlphaFold import AlphaFold
from evodesign.Fitness.Gdt import Gdt
from evodesign.GA.Selection.Tournament import Tournament
from evodesign.GA.Recombination.MiddlePointCrossover import MiddlePointCrossover
from evodesign.GA.Mutation.Switch import Switch
from evodesign.Workspace import Workspace

settings = {
  'Algorithms.PGPD': {
    'numGenerations': 1000,
    'popSize': 100,
    'elitismSize': 20,
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
        'numSelectedCouples': 100,
        'tournamentSize': 5
      }
    },
    'recombination': {
      'GA.Recombination.MiddlePointCrossover': {}
    },
    'mutation': {
      'GA.Mutation.Switch': {}
    }
  }
}
workspace = Workspace('foo', 'bar')
algo = Settings.parse(settings)
test_passed = isinstance(algo, PGPD) and \
              isinstance(algo._predictor, AlphaFold) and \
              isinstance(algo._fitness_fn, Gdt) and \
              isinstance(algo._selection, Tournament) and \
              isinstance(algo._recombination, MiddlePointCrossover) and \
              isinstance(algo._mutation, Switch) and \
              type(algo._predictor._fakemsa_script_path) == str and \
              algo._predictor._fakemsa_script_path == 'fakemsa.py' and \
              type(algo._predictor._alphafold_script_path) == str and \
              algo._predictor._alphafold_script_path == 'run_docker.py' and \
              type(algo._predictor._mgnify_database_path) == str and \
              algo._predictor._mgnify_database_path == '/media/my_usb/mgnify' and \
              type(algo._predictor._data_dir) == str and \
              algo._predictor._data_dir == '/home/com/my_usb' and \
              type(algo._fitness_fn._metrics[1].cutoffs) == list and \
              len(algo._fitness_fn._metrics[1].cutoffs) == 2 and \
              sum([ type(x) == float for x in algo._fitness_fn._metrics[1].cutoffs ]) == 2 and \
              algo._fitness_fn._metrics[1].cutoffs == [ 0.5, 1. ] and \
              type(algo._selection._num_selected_couples) == int and \
              algo._selection._num_selected_couples == 100 and \
              type(algo._selection._tournament_size) == int and \
              algo._selection._tournament_size == 5 and \
              type(algo._num_generations) == int and \
              algo._num_generations == 1000 and \
              type(algo._pop_size) == int and \
              algo._pop_size == 100 and \
              type(algo._elitism_size) == int and \
              algo._elitism_size == 20
print('PASSED' if test_passed else 'FAILED')
