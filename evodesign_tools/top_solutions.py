from argparse import ArgumentParser
import evodesign.Operators as op
import pandas as pd
import os
import json





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

settings_path = f'{args.workspace}/settings.json'
with open(settings_path, 'rt', encoding='utf-8') as json_file:
  settings = json.load(json_file)
k = next(settings.keys())
sort_columns = settings[k]['sortColumns']
sort_ascending = settings[k]['sortAscending']

populations_dir = f'{args.workspace}/populations'
pop_filepaths = os.listdir(populations_dir)
pop_filepaths.sort()
pop = pd.read_csv(f'{populations_dir}/{pop_filepaths[0]}')
pop.sort_values(by=sort_columns, 
                ascending=sort_ascending, 
                inplace=True,
                ignore_index=True)
solutions = pop.iloc[:args.num_solutions].copy()
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
                       args.num_solutions)
solutions.to_csv(f'{args.workspace}/top.csv', index=False)
