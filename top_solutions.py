import json
import os
import pandas as pd
import evodesign.Utils as Utils
from argparse import ArgumentParser





"""
WARNING: this function is NOT working correctly as of yet
"""
def top_solutions_from_all_populations(workspace_dir: str,
                                       num_solutions: int
                                       ) -> pd.DataFrame:
  settings_path = f'{workspace_dir}/settings.json'
  with open(settings_path, 'rt', encoding='utf-8') as json_file:
    settings = json.load(json_file)
  k = list(settings.keys())[0]
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
    solutions = Utils.merge(solutions, 
                            pop, 
                            sort_columns, 
                            sort_ascending, 
                            num_solutions)
  return solutions



if __name__ == '__main__':
  parser = ArgumentParser(prog='top_solutions',
                          description='A utility script for extracting the top '
                          'solutions from all populations produced by an '
                          'evolutionary algorithms')
  parser.add_argument('workspace',
                      help='path to the folder that contains the workspace of '
                      'the evolutionary algorithm which solutions will be '
                      'extracted')
  parser.add_argument('num_solutions',
                      help='the number of solutions that will be extracted')
  args = parser.parse_args()
  solutions = tools.top_solutions_from_all_populations(args.workspace, 
                                                       args.num_solutions)
  solutions.to_csv(f'{args.workspace}/top.csv', index=False)
