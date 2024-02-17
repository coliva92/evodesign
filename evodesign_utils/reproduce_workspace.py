from argparse import ArgumentParser
import pandas as pd
import shutil
import os
import glob





parser = ArgumentParser(prog='reproduce_workspace',
                        description='A utility script for replicating an '
                                    'EvoDesign workspace in order to reproduce '
                                    'the corresponding experiment')
parser.add_argument('source_workspace',
                    help='path to the folder that contains the original '
                         'workspace to be reproduced')
parser.add_argument('destination_workspace',
                    help='path to the folder where the results from the '
                         'reproduced experiment will be stored')
args = parser.parse_args()

pdbs_dir = f'{args.destination_workspace}/pdbs'
os.makedirs(pdbs_dir)
for filepath in glob.glob(f'{args.source_workspace}/pdbs/prot_0001_*.pdb'):
  shutil.copy(filepath, pdbs_dir)

populations_dir = f'{args.destination_workspace}/populations'
os.makedirs(populations_dir)
shutil.copy(f'{args.source_workspace}/populations/pop_0001.csv', 
            populations_dir)
# remove fitness columns
df = pd.read_csv(f'{populations_dir}/pop_0001.csv')
df = df[[ 'generation_id', 'sequence_id', 'sequence', 'survivor' ]]
df.to_csv(f'{populations_dir}/pop_0001.csv')

shutil.copy(f'{args.source_workspace}/initial_rng_state.json',
            args.destination_workspace)

shutil.copy(f'{args.source_workspace}/settings.json', 
            args.destination_workspace)

for filepath in glob.glob(f'{args.source_workspace}/*.pdb'):
  shutil.copy(filepath, args.destination_workspace)
