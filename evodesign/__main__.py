from argparse import ArgumentParser
from .Workspace import Workspace
from .Settings import Settings
from .Exceptions import *
from requests.exceptions import ConnectTimeout
import sys
import json





parser = ArgumentParser(prog='evodesign',
                        description='A rudimentary suite of evolutionary '
                                    'algorithms for protein design.')
parser.add_argument('target_pdb',
                    help='path to the PDB file that contains the target '
                         'protein backbone')
parser.add_argument('settings_filename', 
                    help='path to the JSON file describing the configuration '
                         'of the evolutionary algorithm to be executed')
parser.add_argument('workspace_root',
                    help='path to the folder where all the output files will '
                         'be written')
parser.add_argument('-n', '--num-generations', 
                    type=int, 
                    default=0,
                    help='stops algorithm execution after a the specified '
                         'amount of generations are produced, instead of '
                         'the amount specified in the settings file')
args = parser.parse_args()

with open(args.settings_filename, 'rt', encoding='utf-8') as json_file:
  settings = json.load(json_file)
algorithm = Settings.parse(settings)
reference, population = algorithm.setup(args.target_pdb, args.workspace_root)
Workspace.instance().save_commit_hash()

while True:
  try:
    algorithm(reference, population, numGenerations=args.num_generations)
    break
  except (KeyboardInterrupt, HttpBadRequest, HttpUnknownError):
    print(f'\nINTERRUPTED.\n'
          f'Run `python -m evodesign {args.workspace_root} '
          f'{args.target_pdb} {args.settings_filename}` to resume later.')
    sys.exit(130) # SIGINT
  except (HttpInternalServerError, 
          HttpGatewayTimeout,
          HttpForbidden,
          ConnectTimeout):
    pass
