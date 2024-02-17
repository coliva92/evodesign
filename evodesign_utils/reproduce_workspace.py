from argparse import ArgumentParser
from ..evodesign.Population import Population
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
pop = Population.create_random(1, 4)
df = pd.read_csv(f'{populations_dir}/pop_0001.csv')
df = df[pop.columns]
df.to_csv(f'{populations_dir}/pop_0001.csv')

shutil.copy(f'{args.source_workspace}/initial_rng_state.json',
            args.destination_workspace)
