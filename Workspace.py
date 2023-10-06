from .Individual import Individual
from .Population import Population
from .Statistics import Statistics
import evodesign.FileIO as FileIO
import matplotlib.pyplot as plt
import json
import csv
import os
import random
import time





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
    self.graph_filename = os.path.join(rootFolder, 'fitness_diversity.png')
    self.populations_folder = os.path.join(rootFolder, 'populations')
    self.pdbs_folder = os.path.join(rootFolder, 'pdbs')
    self._seed = time.time()
    random.seed(self._seed)



  def save_algorithm_settings(self, settings: dict) -> None:
    os.makedirs(self.root_folder, exist_ok=True)
    with open(self.settings_filename, 'wt', encoding='utf-8') as json_file:
      json_file.write(json.dumps(settings, indent=2) + '\n')
  


  def save_rng_settings(self) -> None:
    os.makedirs(self.root_folder, exist_ok=True)
    state = random.getstate()
    settings = {
      'seed': self._seed,
      'state': ( state[0], list(state[1]), state[2])
    }
    with open(self.rng_settings_filename, 'wt', encoding='utf-8') as json_file:
      json.dump(settings, json_file)
  


  def load_rng_settings(self) -> None:
    if not os.path.isfile(self.rng_settings_filename):
      self._seed = time.time()
      random.seed(self._seed)
      return
    with open(self.rng_settings_filename, 'rt', encoding='utf-8') as json_file:
      settings = json.load(json_file)
    self._seed = settings['seed']
    random.seed(self._seed)
    if 'state' in settings and settings['state'] is not None:
      settings['state'][1] = tuple(settings['state'][1])
      random.setstate(tuple(settings['state']))



  def save_population(self, population: Population) -> None:
    os.makedirs(self.populations_folder, exist_ok=True)
    filename = population.get_filename(self.populations_folder)
    FileIO.save_population_csv(population, filename)



  def load_latest_population(self) -> Population:
    if not os.path.isdir(self.populations_folder):
      return Population()
    filenames = os.listdir(self.populations_folder)
    if not filenames:
      return Population()
    filename = os.path.join(self.populations_folder, sorted(filenames)[-1])
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
  


  def plot(self) -> None:
    if not os.path.isfile(self.stats_filename):
      return
    data = {
      'iteration_id': [],
      'min_fitness': [],
      'fitness_mean': [],
      'max_fitness': [],
      'sequence_diversity': [],
      'residue_diversity': [],
      'best_sequence_fitness': []
    }
    with open(self.stats_filename, 'rt', encoding='utf-8') as csv_file:
      for row in csv.DictReader(csv_file, dialect='unix'):
        for key in data:
          data[key].append(float(row[key]))
    fig, ax = plt.subplots(ncols=3, figsize=(21, 6))
    fig.suptitle(self.root_folder)
    ax[0].plot(data['iteration_id'], data['fitness_mean'], label='Fitness mean')
    ax[0].fill_between(data['iteration_id'], 
                       data['min_fitness'], 
                       data['max_fitness'], 
                       alpha=0.1)
    ax[0].plot(data['iteration_id'], 
               data['best_sequence_fitness'],
               label='Best solution found')
    ax[0].set_xlabel('Iterations')
    ax[0].set_ylabel('Fitness')
    ax[0].legend(loc='best')
    ax[1].plot(data['iteration_id'], data['sequence_diversity'], color='C2')
    ax[1].set_xlabel('Iterations')
    ax[1].set_ylabel('Population diversity')
    ax[2].plot(data['iteration_id'], data['residue_diversity'], color='C3')
    ax[1].set_xlabel('Iterations')
    ax[1].set_ylabel('Amino acid diversity')
    fig.savefig(self.graph_filename)
