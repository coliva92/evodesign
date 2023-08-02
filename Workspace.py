from typing import List, Optional, Tuple
from .Individual import Individual
import json
import os





class Workspace:
  
  def __init__(self, 
               name: str, 
               targetPdbFilename: str,
               populationFilenames: Optional[List[str]] = None
               ) -> None:
    if populationFilenames is None: populationFilenames = []
    self.reference_filename = targetPdbFilename
    self.name = name
    self.settings_filename = os.path.abspath(os.path.join(self.name, 
                                                       'settings.json'))
    self.stats_filename = self.settings_filename.replace('settings.json', 
                                                      'statistics.csv')
    self.children_filename = self.settings_filename.replace('settings.json', 
                                                         '~children.tmp')
    self.graph_filename = self.settings_filename.replace('settings.json', 
                                                      'fitness.png')
    self.populations_folder = self.settings_filename.replace('settings.json', 
                                                          'populations')
    self.pdbs_folder = self.settings_filename.replace('settings.json', 
                                                   'pdbs')
    self.population_filenames = populationFilenames
    self.memento = {}



  def save_population(self, 
                      iterationId: int, 
                      population: List[Individual]
                      ) -> None:
    os.makedirs(self.populations_folder, exist_ok=True)
    memento = [ individual.as_dict() for individual in population ]
    filename = os.path.join(self.populations_folder, f'pop_{iterationId}.json')
    new_save_file = not os.path.isfile(filename)
    with open(filename, 'wt', encoding='utf-8') as the_file:
      the_file.write(json.dumps(memento, indent=2) + '\n')
    if new_save_file:
      self.population_filenames.append(filename)
      with open(self.settings_filename, 'wt', encoding='utf-8') as the_file:
        the_file.write(json.dumps(self.memento, indent=2) + '\n')



  def load_latest_population(self) ->Tuple[int, List[Individual]]:
    if not self.population_filenames:
      return 0, []
    filename = self.population_filenames[-1]
    with open(filename, 'rt', encoding='utf-8') as json_file:
      individuals_data = json.load(json_file)
    iterationId = len(self.population_filenames) - 1
    population = [ Individual(**params) for params in individuals_data ]
    return iterationId, population



  def backup_children(self, 
                      children: List[Individual]
                      ) -> None:
    memento = [ child.as_dict() for child in children ]
    with open(self.children_filename, 'wt', encoding='utf-8') as json_file:
      json_file.write(json.dumps(memento, indent=2) + '\n')



  def delete_children_backup(self) -> None:
    with open(self.children_filename, 'wt', encoding='utf-8') as the_file:
      the_file.write('[]\n')



  def restore_children_from_backup(self) -> List[Individual]:
    if not os.path.isfile(self.children_filename):
      return []
    with open(self.children_filename, 'rt', encoding='utf-8') as the_file:
      memento = json.load(the_file)
    return [ Individual(**params) for params in memento ]
