from typing import List, Optional
from .Individual import Individual
from .Population import Population
from .Statistics import Statistics
import evodesign.FileIO as FileIO
import matplotlib.pyplot as plt
import json
import csv
import os
import random





class Workspace:
  
  def __init__(self, 
               rootFolder: str,
               targetPdbFilename: str
               ) -> None:
    self.reference_filename = targetPdbFilename
    self.root_folder = rootFolder
    self.settings_filename = os.path.join(rootFolder, 'settings.json')
    self.rng_settings_filename = os.path.join(rootFolder, 'rng.json')
    self.stats_filename = os.path.join(rootFolder, 'statistics.csv')
    self.children_filename = os.path.join(rootFolder, '~children.tmp')
    self.graphs_filename = os.path.join(rootFolder, 'fitness_diversity.png')
    self.populations_folder = os.path.join(rootFolder, 'populations')
    self.pdbs_folder = os.path.join(rootFolder, 'pdbs')
    self.algorithm_settings = None



  def save_algorithm_settings(self) -> None:
    with open(self.settings_filename, 'wt', encoding='utf-8') as json_file:
      json_file.write(json.dumps(self.algorithm_settings, indent=2) + '\n')
  


  def save_rng_settings(self, state: tuple) -> None:
    state = ( state[0], list(state[1]), state[2])
    with open(self.rng_settings_filename, 'wt', encoding='utf-8') as json_file:
      json.dump(state, json_file)
  


  def load_rng_settings(self) -> bool:
    if not os.path.isfile(self.rng_settings_filename):
      return False
    with open(self.rng_settings_filename, 'rt', encoding='utf-8') as json_file:
      state = json.load(json_file)
    state[1] = tuple(state[1])
    random.setstate(tuple(state))
    return True



  def save_population(self, population: Population) -> None:
    os.makedirs(self.populations_folder, exist_ok=True)
    filename = population.get_filename(self.populations_folder)
    FileIO.save_population_csv(population, filename)



  def load_latest_population(self) -> Population:
    filenames = os.listdir(self.populations_folder)
    if not filenames:
      return Population()
    filename = sorted(filenames)[-1]
    return FileIO.load_population_csv(filename, len(filenames))



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
  


  def plot_fitness(self) -> None:
    if not os.path.isfile(self.stats_filename):
      return
    data = {
      'iteration_id': [],
      'min_fitness': [],
      'fitness_mean': [],
      'max_fitness': [],
      'best_sequence_fitness': []
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
    ax.set_title(self.root_folder)
    ax.set_xlabel('Iterations')
    ax.set_ylabel('Fitness')
    ax.legend(loc='best')
    fig.savefig(self.graph_filename)
