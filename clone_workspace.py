import os
import shutil
import glob
import pandas as pd
from argparse import ArgumentParser





def reproduce_workspace(source_dir: str, 
                        destination_dir: str
                        ) -> None:
  pdbs_dir = f'{destination_dir}/pdbs'
  os.makedirs(pdbs_dir)
  for pdb_path in glob.glob(f'{source_dir}/pdbs/prot_0001_*.pdb'):
    shutil.copy(pdb_path, pdbs_dir)
  populations_dir = f'{destination_dir}/populations'
  os.makedirs(populations_dir)
  shutil.copy(f'{source_dir}/populations/pop_0001.csv', 
              populations_dir)
  # remove fitness columns
  df = pd.read_csv(f'{populations_dir}/pop_0001.csv')
  df = df[[ 'generation_id', 'sequence_id', 'sequence', 'survivor' ]]
  df.to_csv(f'{populations_dir}/pop_0001.csv', index=False)
  shutil.copy(f'{source_dir}/initial_rng_state.json', destination_dir)
  shutil.copy(f'{source_dir}/settings.json', destination_dir)
  for pdb_path in glob.glob(f'{source_dir}/*.pdb'):
    shutil.copy(pdb_path, destination_dir)



if __name__ == '__main__':
  parser = ArgumentParser(prog='clone_workspace',
                          description='A utility script for replicating an '
                                      'EvoDesign workspace in order to '
                                      'reproduce its execution')
  parser.add_argument('source_workspace',
                      help='path of the original workspace to be reproduced')
  parser.add_argument('destination_workspace',
                      help='path to the folder where the workspace clone will be stored')
  args = parser.parse_args()
  reproduce_workspace(args.source_workspace, args.destination_workspace)
