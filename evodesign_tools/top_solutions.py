from argparse import ArgumentParser
import evodesign_tools as tools





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
