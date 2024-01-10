from argparse import ArgumentParser
import evodesign.Settings as Settings
from .Exceptions import *
from requests.exceptions import ConnectTimeout
import sys





parser = ArgumentParser(prog='evodesign',
                        description='A rudimentary framework for testing ' +
                                    'evolutionary algorithms for protein ' + 
                                    'design.')
parser.add_argument('settings_filename', 
                    help='name or path for the JSON configuration file for '+ 
                    'the evolutionary algorithm to be ran.')
parser.add_argument('-a', '--add-iterations', type=int, default=0,
                    help='runs the algorithm for a number of iterations that ' +
                    'is added to the value specified in the `numIterations` ' +
                    'field in the specified settings file.')
args = parser.parse_args()
filename = args.settings_filename
while True:
  try:
    algorithm = Settings.load_algorithm_from_settings(filename, 
                                                      args.add_iterations)
    if filename != algorithm.workspace.settings_filename:
      algorithm.workspace.save_algorithm_settings(algorithm.as_json())
      filename = algorithm.workspace.settings_filename
    algorithm.workspace.load_rng_json()
    population = algorithm.workspace.load_latest_population()
    algorithm(population)
    algorithm.workspace.plot()
    print(f'\nCOMPLETED.\n' +
          f'Best sequence found: {algorithm.best_solution.sequence}\n' + 
          f'Fitness: {algorithm.best_solution.fitness:0.4f}')
    break
  except (KeyboardInterrupt, HttpUnknownError):
    print(f'\nINTERRUPTED.\n' +
          f'Run `python -m evodesign {filename}` to resume later.')
    algorithm.workspace.plot()
    sys.exit(130)
  except (HttpInternalServerError, 
          HttpGatewayTimeout,
          HttpForbidden,
          ConnectTimeout
          ) as e:
    pass
