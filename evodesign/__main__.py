from argparse import ArgumentParser
import evodesign.Settings as Settings
from .Exceptions import *
from requests.exceptions import ConnectTimeout
import sys
import json
from .Context import Context





parser = ArgumentParser(prog='evodesign',
                        description='A rudimentary suite of evolutionary '
                                    'algorithms for protein design.')
parser.add_argument('target_pdb_path',
                    help='path to the PDB file that contains the target '
                         'protein backbone')
parser.add_argument('settings_path', 
                    help='path to the JSON file describing the configuration '
                         'of the evolutionary algorithm to be executed')
parser.add_argument('workspace_root',
                    help='path to the folder where all the output files will '
                         'be written')
parser.add_argument('-n', '--num_generations', 
                    type=int, 
                    default=0,
                    help='stops algorithm execution after a the specified '
                         'amount of generations are produced, instead of '
                         'the amount specified in the settings file')
parser.add_argument('-f', '--target_fasta_path',
                    type=str,
                    default=None,
                    help='path to the FASTA file containing the target '
                         'amino acid sequence for when the goal is to find '
                         'a different sequence with the same properties as the '
                         'target protein')
parser.add_argument('-r', '--sequence_restrictions',
                    type=str,
                    default=None,
                    help='path to the JSON file describing the allowed amino acids '
                         'for certain positions in the designed sequences')
args = parser.parse_args()

context = Context.create(args.target_pdb_path, 
                         args.target_fasta_path, 
                         args.num_generations,
                         args.sequence_restrictions)
with open(args.settings_path, 'rt', encoding='utf-8') as json_file:
    settings = json.load(json_file)
algorithm = Settings.parse(settings)
population = algorithm.setup(context, args.workspace_root)

while True:
    try:
        algorithm(population)
        break
    except (KeyboardInterrupt, HttpBadRequest, HttpUnknownError):
        temp = f'-s {args.target_fasta_path} ' if args.target_fasta_path else ''
        print(f'\nINTERRUPTED.\n'
              f'Run `python -m evodesign {args.workspace_root} '
              f'{args.target_pdb_path} {args.settings_path}` '
              f'{temp} to resume later.')
        sys.exit(130) # SIGINT
    except (HttpInternalServerError, 
            HttpGatewayTimeout,
            HttpForbidden,
            ConnectTimeout):
        pass
