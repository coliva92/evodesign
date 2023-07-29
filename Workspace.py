from typing import List
from .Population import Individual
import json
import os





class Workspace:
  """La interfaz para manejar el acceso a la carpeta donde se almacenan los 
  datos producidos por el algoritmo evolutivo.
  """
  
  def __init__(self, 
               name: str, 
               targetPdbFilename: str
               ) -> None:
    """Constructor.
    - `name`: el nombre de la carpeta.
    - `targetPdbFilename`: el nombre del archivo PDB que contiene la estructura 
      objetivo.
    """
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
    self.population_filenames = []
    self.memento = {}



  def save_population(self, 
                      iterationId: int, 
                      population: List[Individual]
                      ) -> None:
    """Respalda el conjunto de individuos especificado por `population` en un 
    archivo JSON localizado en `{self.name}/populations/pop_{iterationId}.json`.
    - `iterationId`: el identificador de la iteración que corresponde a la 
    población a respaldar.
    - `population`: la colección de individuos que van a respaldarse. 
    """
    os.makedirs(self.populations_folder, exist_ok=True)
    memento = [ individual.get_memento() for individual in population ]
    filename = os.path.join(self.populations_folder, f'pop_{iterationId}.json')
    new_save_file = not os.path.isfile(filename)
    with open(filename, 'wt', encoding='utf-8') as the_file:
      the_file.write(json.dumps(memento, indent=2) + '\n')
    if new_save_file:
      self.population_filenames.append(filename)
      with open(self.settings_filename, 'wt', encoding='utf-8') as the_file:
        the_file.write(json.dumps(self.memento, indent=2) + '\n')



  def backup_children(self, 
                      children: List[Individual]
                      ) -> None:
    """Respalda la lista de hijos especificada por `children` en un archivo 
    JSON localizado en `{self.name}/~children.tmp`.
    - `children`: la lista de secuencias hijas que van a respaldarse.
    """
    memento = [ individual.get_memento() for individual in children ]
    with open(self.children_filename, 'wt', encoding='utf-8') as the_file:
      the_file.write(json.dumps(memento, indent=2) + '\n')
  


  def update_children_backup(self, 
                             children: List[Individual]
                             ) -> None:
    memento = [ child.get_memento() for child in children ]
    with open(self.children_filename, 'wt', encoding='utf-8') as the_file:
      the_file.write(json.dumps(memento, indent=2) + '\n')



  def delete_children_backup(self) -> None:
    """Vacía el contenido del archivo `{self.name}/~children.tmp`, usado para 
    respaldar los datos de las secuencias hijas producidas.
    """
    with open(self.children_filename, 'wt', encoding='utf-8') as the_file:
      the_file.write('[]\n')



  def restore_children_from_backup(self) -> List[Individual]:
    """Carga a un arreglo el contenido del archivo `{self.name}/~children.tmp`, 
    usado para respaldar los datos de las secuencias hijas producidas.
    """
    if not os.path.isfile(self.children_filename):
      return []
    with open(self.children_filename, 'rt', encoding='utf-8') as the_file:
      memento = json.load(the_file)
    return [ Individual(**params) for params in memento ]
