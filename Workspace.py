from typing import List, Optional, Callable
from .Individual import Individual
from .Population import Population
from .Statistics import Statistics
import evodesign.FileIO as FileIO
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
    self.stats_filename = os.path.join(self.name, 'statistics.csv')
    self.children_filename = os.path.join(self.name, '~children.tmp')
    self.graph_filename = os.path.join(self.name, 'fitness.png')
    self.populations_folder = os.path.join(self.name, 'populations')
    self.pdbs_folder = self.settings_filename.replace('settings.json', 'pdbs')
    self.population_filenames = populationFilenames
    self.json_factory = jsonFactory



  def save_population_and_update_settings(self, population: Population) -> None:
    os.makedirs(self.populations_folder, exist_ok=True)
    filename = population.get_filename(self.populations_folder)
    if FileIO.save_population_csv(population, filename):
      settings_json = self.json_factory()
      self.population_filenames.append(filename)
      settings_json['__savedPopulations'] = self.population_filenames
      with open(self.settings_filename, 'wt', encoding='utf-8') as json_file:
        json_file.write(json.dumps(settings_json, indent=2) + '\n')



  def load_latest_population(self) -> Population:
    if not self.population_filenames:
      return Population()
    filename = self.population_filenames[-1]
    return FileIO.load_population_csv(filename, len(self.population_filenames))



  def backup_children(self, children: Population) -> None:
    FileIO.save_population_json(children, self.children_filename)



  def restore_children_from_backup(self) -> Population:
    if not os.path.isfile(self.children_filename):
      return Population()
    return FileIO.load_population_json(self.children_filename)
  


  def delete_children_backup(self) -> None:
    with open(self.children_filename, 'wt', encoding='utf-8') as json_file:
      json_file.write('[]\n')
  


  def save_statistics(self, 
                      stats: Statistics, 
                      best_solution: Individual) -> None:
    file_exists = os.path.isfile(self.stats_filename)
    data = stats.as_dict()
    data['best_sequence_fitness'] = best_solution.fitness
    data['best_sequence'] = best_solution.sequence
    with open(self.stats_filename, 'at', encoding='utf-8') as csv_file:
      writer = csv.DictWriter(csv_file, 
                              data.keys(), 
                              dialect='unix',
                              quoting=csv.QUOTE_NONE)
      if not file_exists: 
        writer.writeheader()
      writer.writerow(data)
    return
  


  def plot_fitness(self) -> None:
    if not os.path.isfile(self.stats_filename):
      return
    data = {
      'iteration_id': [],
      'min_fitness': [],
      'fitness_mean': [],
      'max_fitness': []
    }
    with open(self.stats_filename, 'rt', encoding='utf-8') as csv_file:
      for row in csv.DictReader(csv_file, dialect='unix'):
        for key in data:
          data[key].append(float(row[key]))
    fig, ax = plt.subplots(1)
    ax.plot(data['iteration_id'], data['fitness_mean'], label='Fitness mean')
    ax.fill_between(data['iteration_id'], 
                    data['min_fitness'], 
                    data['max_fitness'], 
                    alpha=0.1)
    ax.plot(data['iteration_id'], 
            data['best_sequence_fitness'], 
            label='Best solution found')
    ax.set_title(self.name)
    ax.set_xlabel('Iterations')
    ax.set_ylabel('Fitness')
    ax.legend(loc='best')
    fig.savefig(self.graph_filename)
