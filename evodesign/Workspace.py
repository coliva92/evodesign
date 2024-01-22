from Individual import Individual
from Population import Population
from Statistics import Statistics
import FileIO
import matplotlib.pyplot as plt
import json
import csv
import os
import random
import time





class Workspace:

  _instance = None



  @classmethod
  def instance(self):
    """
    Returns
    -------
    Workspace
        The singleton instance for the Workspace object.
    """
    return self._instance



  def __new__(cls, *args, **kwargs):
    if not cls._instance:
      cls._instance = super(Workspace, cls).__new__(cls)
    return cls._instance



  def __init__(self, 
               root: str,
               targetPdb: str
               ) -> None:
    """
    Interface for managing the file system folder where all the output files 
    will be stored.

    Parameters
    ----------
    root : str
        The folder where all the output files and folders will be stored.
    targetPdb : str
        The path to the PDB file containing the backbone for which an 
        amino acid sequence will be designed.
    """
    self.target_pdb = targetPdb
    self.root = root
    self.settings = f'{self.root}/settings.json'
    self.rng_initial_state = f'{self.root}/rng_initial_state.json'
    self.rng_checkpoint = f'{self.root}/.rng_checkpoint'
    self.statistics = f'{self.root}/statistics.csv'
    self.generation_checkpoint = f'{self.root}/.generation_checkpoint'
    self.fitness_diversity = f'{self.root}/fitness_diversity.png'
    self.populations = f'{self.root}/populations'
    self.pdbs = f'{self.root}/pdbs'
    # ----
    self._seed = time.time()
    random.seed(self._seed)



  def save_algorithm_settings(self, settings: dict) -> None:
    os.makedirs(self.root_folder, exist_ok=True)
    with open(self.settings_filename, 'wt', encoding='utf-8') as json_file:
      json_file.write(json.dumps(settings, indent=2) + '\n')
  


  def save_rng_json(self, checkpoint: bool = True) -> None:
    os.makedirs(self.root_folder, exist_ok=True)
    filename = self.rng_checkpoint_filename \
               if checkpoint \
               else self.rng_settings_filename
    state = random.getstate()
    settings = {
      'seed': self._seed,
      'state': ( state[0], list(state[1]), state[2] )
    }
    with open(filename, 'wt', encoding='utf-8') as json_file:
      json.dump(settings, json_file)
  


  def load_rng_json(self, checkpoint: bool = True) -> None:
    filename = self.rng_checkpoint_filename \
               if checkpoint \
               else self.rng_settings_filename
    if not os.path.isfile(filename):
      self._seed = time.time()
      random.seed(self._seed)
      self.save_rng_json(checkpoint=False)
      return
    with open(filename, 'rt', encoding='utf-8') as json_file:
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
    filenames = sorted(os.listdir(self.populations_folder))
    if not filenames:
      return Population()
    filename = os.path.join(self.populations_folder, filenames[-1])
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
      'sequence_identity': [],
      'residue_identity': [],
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
    ax[1].plot(data['iteration_id'], data['sequence_identity'], color='C2')
    ax[1].set_xlabel('Iterations')
    ax[1].set_ylabel('Population diversity')
    ax[2].plot(data['iteration_id'], data['residue_identity'], color='C3')
    ax[2].set_xlabel('Iterations')
    ax[2].set_ylabel('Amino acid diversity')
    fig.savefig(self.graph_filename)
