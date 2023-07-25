from argparse import ArgumentParser
import evodesign.Statistics as Statistics





parser = ArgumentParser()
parser.add_argument('filename')
parser.add_argument('figFilename')
parser.add_argument('figTitle')
args = parser.parse_args()
Statistics.plot_fitness_over_iterations(args.filename, 
                                        args.figFilename, 
                                        args.fitTitle)
