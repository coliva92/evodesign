import evodesign.Operators as op
import pandas as pd
import glob
import shutil
import json
import os





def reproduce_workspace(source_dir: str, 
                        destination_dir: str
                        ) -> None:
  pdbs_dir = f'{destination_dir}/pdbs'
  os.makedirs(pdbs_dir)
  for filepath in glob.glob(f'{source_dir}/pdbs/prot_0001_*.pdb'):
    shutil.copy(filepath, pdbs_dir)
  populations_dir = f'{destination_dir}/populations'
  os.makedirs(populations_dir)
  shutil.copy(f'{source_dir}/populations/pop_0001.csv', 
              populations_dir)
  # remove fitness columns
  df = pd.read_csv(f'{populations_dir}/pop_0001.csv')
  df = df[[ 'generation_id', 'sequence_id', 'sequence', 'survivor' ]]
  df.to_csv(f'{populations_dir}/pop_0001.csv')
  shutil.copy(f'{source_dir}/initial_rng_state.json', destination_dir)
  shutil.copy(f'{source_dir}/settings.json', destination_dir)
  for filepath in glob.glob(f'{source_dir}/*.pdb'):
    shutil.copy(filepath, destination_dir)



def top_solutions_from_all_populations(workspace_dir: str,
                                       num_solutions: int
                                       ) -> pd.DataFrame:
  settings_path = f'{workspace_dir}/settings.json'
  with open(settings_path, 'rt', encoding='utf-8') as json_file:
    settings = json.load(json_file)
  k = next(settings.keys())
  sort_columns = settings[k]['sortColumns']
  sort_ascending = settings[k]['sortAscending']
  populations_dir = f'{workspace_dir}/populations'
  pop_filepaths = os.listdir(populations_dir)
  pop_filepaths.sort()
  pop = pd.read_csv(f'{populations_dir}/{pop_filepaths[0]}')
  pop.sort_values(by=sort_columns, 
                  ascending=sort_ascending, 
                  inplace=True,
                  ignore_index=True)
  solutions = pop.iloc[:num_solutions].copy()
  for filepath in pop_filepaths[1:]:
    pop = pd.read_csv(f'{populations_dir}/{filepath}')
    pop.sort_values(by=sort_columns, 
                    ascending=sort_ascending, 
                    inplace=True,
                    ignore_index=True)
    solutions = op.merge(solutions, 
                        pop, 
                        sort_columns, 
                        sort_ascending, 
                        num_solutions)
  return solutions
