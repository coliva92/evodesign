from argparse import ArgumentParser
import evodesign.Settings as Settings
from .Exceptions import *
import sys
from requests.exceptions import ConnectTimeout





parser = ArgumentParser(prog='evodesign',
                        description='A rudimentary framework for testing ' +
                                    'evolutionary algorithms for protein ' + 
                                    'design.')
parser.add_argument('settings_filename', 
                    help='name or path for the JSON configuration file for '+ 
                    'the evolutionary algorithm to be ran.')
args = parser.parse_args()
filename = args.settings_filename
while True:
  try:
    algorithm = Settings.load_algorithm_from_settings(filename)
    filename = algorithm.workspace.settings_filename
    population = algorithm.workspace.load_latest_population()
    algorithm(population)
    algorithm.workspace.plot_fitness()
    print(f'COMPLETED.\n' +
          f'Best sequence found: {algorithm.best_solution.sequence}\n' + 
          f'Fitness: {algorithm.best_solution.fitness:0.4f}')
    break
  except (KeyboardInterrupt, HttpUnknownError):
    print(f'\nINTERRUPTED.\n' +
          f'Run `python -m evodesign {filename}` to resume later.')
    algorithm.workspace.plot_fitness()
    sys.exit(130)
  except (HttpInternalServerError, 
          HttpGatewayTimeout, 
          HttpForbidden, 
          ConnectTimeout
          ) as e:
    pass
