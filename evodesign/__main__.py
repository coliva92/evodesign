from argparse import ArgumentParser
import Settings
from Exceptions import *
from requests.exceptions import ConnectTimeout
import sys





parser = ArgumentParser(prog='evodesign',
                        description='A rudimentary suite of evolutionary ' + 
                                    'algorithms for protein design.')
parser.add_argument('target_pdb',
                    help='path to the PDB file that contains the target ' + 
                         'protein backbone')
parser.add_argument('settings_filename', 
                    help='path to the JSON file describing the configuration ' +
                         'of the evolutionary algorithm to be executed')
parser.add_argument('workspace_root',
                    help='path to the folder where all the output files will ' +
                         'be written')
parser.add_argument('-n', '--for-generations', 
                    type=int, 
                    default=0,
                    help='stops algorithm execution after a the specified ' + 
                         'amount of generations are produced, instead of ' +
                         'the amount specified by the `generations` field ' + 
                         'in the settings file')
args = parser.parse_args()


while True:
  try:
    algorithm = Settings.load_algorithm_from_settings(filename, 
                                                      args.add_iterations)
    if filename != algorithm.workspace.settings_filename:
      algorithm.workspace.save_algorithm_settings(algorithm.settings())
      filename = algorithm.workspace.settings_filename
    algorithm.workspace.load_rng_json()
    population = algorithm.workspace.load_latest_population()
    algorithm(population)
    algorithm.workspace.plot()
    break


  except (KeyboardInterrupt, HttpUnknownError):
    print(f'\nINTERRUPTED.\n' +
          f'Run `python -m evodesign {args.workspace_root} ' +
          f'{args.target_pdb} {args.settings_filename}` to resume later.')
    algorithm.workspace.plot()
    sys.exit(130) # SIGINT
  

  except (HttpInternalServerError, 
          HttpGatewayTimeout,
          HttpForbidden,
          ConnectTimeout):
    pass
