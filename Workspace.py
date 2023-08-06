from typing import List, Optional, Callable
from .Individual import Individual
from .Population import Population
from .Statistics import Statistics
import matplotlib.pyplot as plt
import json
import csv
import os





class Workspace:
  
  def __init__(self, 
               name: str, 
               jsonFactory: Callable[[], dict],
               targetPdbFilename: str,
               populationFilenames: Optional[List[str]] = None
               ) -> None:
    if populationFilenames is None: populationFilenames = []
    self.reference_filename = targetPdbFilename
    self.name = name
    self.settings_filename = os.path.join(self.name, 'settings.json')
    self.stats_filename = self.settings_filename.replace('settings.json', 
                                                         'statistics.csv')
    self.children_filename = self.settings_filename.replace('settings.json', 
                                                            '~children.tmp')
    self.graph_filename = self.settings_filename.replace('settings.json', 
                                                         'fitness.png')
    self.populations_folder = self.settings_filename.replace('settings.json', 
                                                             'populations')
    self.pdbs_folder = self.settings_filename.replace('settings.json', 'pdbs')
    self.population_filenames = populationFilenames
    self.json_factory = jsonFactory
  


  def save_population(self, population: Population) -> None:
    os.makedirs(self.populations_folder, exist_ok=True)
    filename = population.get_json_filename(self.populations_folder)
    is_new_file = not os.path.isfile(filename)
    with open(filename, 'wt', encoding='utf-8') as the_file:
      the_file.write(json.dumps(population.as_json(), indent=2) + '\n')
    if is_new_file:
      settings_json = self.json_factory()
      self.population_filenames.append(filename)
      settings_json['__savedPopulations'] = self.population_filenames
      with open(self.settings_filename, 'wt', encoding='utf-8') as the_file:
        the_file.write(json.dumps(settings_json, indent=2) + '\n')



  def load_latest_population(self) -> Population:
    if not self.population_filenames:
      return Population()
    filename = self.population_filenames[-1]
    with open(filename, 'rt', encoding='utf-8') as json_file:
      pop_data = json.load(json_file)
    return Population([ Individual(**params) for params in pop_data ],
                      len(self.population_filenames))



  def backup_children(self, children: Population) -> None:
    with open(self.children_filename, 'wt', encoding='utf-8') as json_file:
      json_file.write(json.dumps(children.as_json(), indent=2) + '\n')



  def restore_children_from_backup(self) -> Population:
    if not os.path.isfile(self.children_filename):
      return Population()
    with open(self.children_filename, 'rt', encoding='utf-8') as json_file:
      backup = json.load(json_file)
    return Population([ Individual(**params) for params in backup ])
  


  def delete_children_backup(self) -> None:
    with open(self.children_filename, 'wt', encoding='utf-8') as json_file:
      json_file.write('[]\n')
  


  def save_statistics(self, 
                      stats: Statistics, 
                      best_solution: Individual) -> None:
    file_exists = os.path.isfile(self.stats_filename)
    data = stats.as_dict()
    data['best_sequence_fitness'] = best_solution.fitness
    data['best_sequence'] = f'"{best_solution.sequence}"'
    with open(self.stats_filename, 'at', encoding='utf-8') as csv_file:
      writer = csv.DictWriter(csv_file, 
                              data.keys(), 
                              dialect='unix', 
                              quotechar="'",
                              quoting=csv.QUOTE_NONE,
                              escapechar='|')
      if not file_exists: 
        writer.writeheader()
      writer.writerow(data)
    return
  


  def plot_fitness(self) -> None:
    if not os.path.isfile(self.stats_filename):
      return
    iterations, minimums, averages, maximums, best_solutions = [],[],[],[],[]
    with open(self.stats_filename, 'rt', encoding='utf-8') as csv_file:
      for row in csv.DictReader(csv_file, dialect='unix'):
        iterations.append(int(row['iteration_id']))
        minimums.append(float(row['min_fitness']))
        averages.append(float(row['fitness_mean']))
        maximums.append(float(row['max_fitness']))
        best_solutions.append(float(row['best_sequence_fitness']))
    fig, ax = plt.subplots(1)
    ax.plot(iterations, averages, label='Fitness mean')
    ax.fill_between(iterations, minimums, maximums, alpha=0.3)
    ax.plot(iterations, best_solutions, label='Best solution found')
    ax.set_title(self.name)
    ax.set_xlabel('Iterations')
    ax.set_ylabel('Fitness')
    ax.legend(loc='best')
    ax.grid()
    fig.savefig(self.graph_filename)
