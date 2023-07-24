from argparse import ArgumentParser
import evodesign.Storage as Storage
import evodesign.Prediction.ESMFold as ESMFold





parser = ArgumentParser(prog='evodesign',
                        description='A rudimentary framework for testing ' +
                                    'evolutionary algorithms for protein ' + 
                                    'design.')
parser.add_argument('setup_filename', 
                    help='name or path for the JSON configuration file for '+ 
                    'the evolutionary algorithm to be ran.')
args = parser.parse_args()
filename = args.setup_filename
while True:
  try:
    memento = Storage.load_memento_from_json_file(filename)
    algorithm = Storage.load_algorithm_from_memento(memento)
    iterationId, population = Storage.load_population_from_memento(memento)
    algorithm.run(iterationId, population)
    break
  except RuntimeError as e:
    ESMFold.handle_api_errors(e)
    filename = algorithm.workspace.setup_filename
    continue
