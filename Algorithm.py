from abc import ABC, abstractmethod
from typing import List, Optional
from .Population import Individual
import evodesign.Chain as Chain
from .Prediction import Predictor
import json
import os





class _Workspace:
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
    self.setup_filename = os.path.abspath(os.path.join(self.name, 
                                                       'settings.json'))
    self.stats_filename = self.setup_filename.replace('settings.json', 
                                                      'statistics.csv')
    self.children_filename = self.setup_filename.replace('settings.json', 
                                                         '~children.tmp')
    self.graph_filename = self.setup_filename.replace('settings.json', 
                                                      'fitness.png')
    self.populations_folder = self.setup_filename.replace('settings.json', 
                                                          'populations')
    self.pdbs_folder = self.setup_filename.replace('settings.json', 
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
      with open(self.setup_filename, 'wt', encoding='utf-8') as the_file:
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





class Algorithm(ABC):
  """La representación genérica de algún algoritmo evolutivo. 
  """
  
  def __init__(self,
               workspaceName: str,
               targetPdbFilename: str,
               predictor: Predictor
               ) -> None:
    """
    Constructor.
    - `workspaceName`: el nombre de la carpeta donde se guardarán los datos 
      producidos durante la ejecución del algoritmo.
    - `targetPdbFilename`: el nombre del archivo PDB que contiene la estructura 
      de la proteína objetivo; esta estructura se utilizará para evaluar la 
      calidad de la proteína diseñada.
    - `predictor`: el algoritmo de predicción de la estructura de una proteína 
      a utilizar; este algoritmo se utiliza principalemnte para corroborar que 
      la secuencia diseñada se pliegue en la estructura objetivo, especificada 
      por `targetPdbFilename`.
    """
    super().__init__()
    self._predictor = predictor
    reference = Chain.load_structure_from_pdb_file(targetPdbFilename)
    self._sequence_length = Chain.count_residues_from_chain(reference)
    self._reference_backbone = Chain.filter_backbone_atoms_from_chain(reference)
    self.workspace = _Workspace(workspaceName, targetPdbFilename)
    self.best_solution = None



  @classmethod
  @abstractmethod
  def get_type(cls) -> str:
    """Retorna el nombre del algoritmo.
    """
    raise NotImplementedError



  @abstractmethod
  def run(self, 
          iterationId: int, 
          population: List[Individual]
          ) -> None:
    """Inicia la ejecución del algoritmo evolutivo.
    - `iterationId`: un entero que indica a partir de qué valor comienzan a 
      contarse las iteraciones (o generaciones) del algoritmo.
    - `population`: la colección de individuos que conforman la población sobre 
      la que comenzará a trabajar el algoritmo.
    """
    raise NotImplementedError



  def _get_params_memento(self) -> dict:
    """Retorna un diccionario que contiene los parámetros de los componentes y 
    que conforman el algoritmo; este diccionario normalmente se utiliza para 
    respaldar en un archivo la configuración del algoritmo, para así poder 
    restaurarla en caso de que se interrumpa la ejecución.
    """
    return {
      'workspaceName': self.workspace.name,
      'targetPdbFilename': self.workspace.reference_filename,
      'predictor': self._predictor.get_name()
    }



  def create_memento(self) -> dict:
    return {
      'algorithmType': self.get_type(),
      'algorithmParams': self._get_params_memento(),
      '__savedPopulations': self.workspace.population_filenames
    }
  


  def append_and_cache_memento(self, 
                               population_filenames: Optional[List[str]] = None
                               ) -> None:
    if population_filenames == None: population_filenames = []
    self.workspace.population_filenames = population_filenames
    self.workspace.memento = self.create_memento()



  def _get_pdb_filename(self, individual: Individual) -> str:
    return os.path.join(self.workspace.pdbs_folder, f'prot_{individual.id}.pdb')
