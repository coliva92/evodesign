from argparse import ArgumentParser
import evodesign_tools as tools





if __name__ == '__main__':
  parser = ArgumentParser(prog='reproduce_workspace',
                          description='A utility script for replicating an '
                                      'EvoDesign workspace in order to '
                                      'reproduce the corresponding experiment')
  parser.add_argument('source_workspace',
                      help='path to the folder that contains the original '
                          'workspace to be reproduced')
  parser.add_argument('destination_workspace',
                      help='path to the folder where the results from the '
                          'reproduced experiment will be stored')
  args = parser.parse_args()
  tools.reproduce_workspace(args.source_workspace, args.destination_workspace)
