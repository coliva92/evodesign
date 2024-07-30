from argparse import ArgumentParser
import evodesign.Utils as Utils





if __name__ == '__main__':
    parser = ArgumentParser(prog='clone_workspace',
                            description='A utility script for replicating an '
                                        'EvoDesign workspace in order to '
                                        'reproduce its execution')
    parser.add_argument('source_workspace',
                        help='path of the original workspace to be reproduced')
    parser.add_argument('destination_workspace',
                        help='path to the folder where the workspace clone '
                             'will be stored')
    args = parser.parse_args()
    Utils.clone_workspace(args.source_workspace, args.destination_workspace)
