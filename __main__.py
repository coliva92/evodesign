from argparse import ArgumentParser
import evodesign.JSON as JSON





parser = ArgumentParser(prog='evodesign',
                        description='A rudimentary framework for testing ' +
                                    'evolutionary algorithms for protein ' + 
                                    'design.')
parser.add_argument('settings_filename', 
                    help='name or path for the JSON configuration file for '+ 
                    'the evolutionary algorithm to be ran.')
args = parser.parse_args()
filename = args.settings_filename
while True:
  try:
    settings = JSON.load_dict_from_json(filename)
    algorithm = JSON.create_algorithm_from_settings_dict(settings)
    iterationId, population = JSON.load_latest_population_from_settings_dict(
      settings)
    algorithm.run(iterationId, population)
    break
  except RuntimeError as e:
    filename = algorithm.workspace.settings_filename
    if e == 'ESMFold API max requests reached' or e == 'Forbidden': 
      print(f'Interrupted.\n' +
            f'Run `python -m evodesign {filename}` to resume later.')
      break
    continue
print(f'Completed.\n' +
      f'Best sequence: {algorithm.best_solution.sequence}\n' + 
      f'Fitness: {algorithm.best_solution.fitness:0.4f}')
