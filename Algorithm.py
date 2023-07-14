from abc import ABC, abstractmethod
from typing import List, Optional
from .Population import Individual
import evodesign.PDB as PDB
import evodesign.Chain as Chain
from .Prediction import Predictor
import json
import os





class _Workspace:
  
  def __init__(self, 
               name: str, 
               referencePdbFilename: str) -> None:
    self.reference_filename = referencePdbFilename
    self.name = name
    self.setup_filename = os.path.abspath(os.path.join(self.name, 
                                                       'setup.json'))
    self.stats_filename = self.setup_filename.replace('setup.json', 
                                                      'statistics.csv')
    self.children_filename = self.setup_filename.replace('setup.json', 
                                                         '~children.json')
    self.graph_filename = self.setup_filename.replace('setup.json', 
                                                      'fitness.png')
    self.populations_folder = self.setup_filename.replace('setup.json', 
                                                          'populations')
    self.pdbs_folder = self.setup_filename.replace('setup.json', 
                                                   'pdbs')
    self.saved_populations = []
    self.memento = {}
    self.children_memento = []



  def save_population(self, 
                      iterationId: int, 
                      population: List[Individual]) -> None:
    os.makedirs(self.populations_folder, exist_ok=True)
    memento = [ individual.get_memento() for individual in population ]
    filename = os.path.join(self.populations_folder, f'pop_{iterationId}.json')
    new_save_file = not os.path.isfile(filename)
    with open(filename, 'wt', encoding='utf-8') as file:
      file.write(json.dumps(memento, indent=2) + '\n')
    if new_save_file:
      self.saved_populations.append(filename)
      with open(self.setup_filename, 'wt', encoding='utf-8') as file:
        file.write(json.dumps(self.memento, indent=2) + '\n')



  def backup_children(self, children: List[Individual]) -> None:
    memento = [ individual.get_memento() for individual in children ]
    with open(self.children_filename, 'wt', encoding='utf-8') as file:
      file.write(json.dumps(memento, indent=2) + '\n')
    self.children_memento = memento
  


  def update_children_backup(self, idx: int, child: Individual) -> None:
    memento = child.get_memento()
    self.children_memento[idx] = memento
    with open(self.children_filename, 'wt', encoding='utf-8') as file:
      file.write(json.dumps(self.children_memento, indent=2) + '\n')
  


  def delete_children_backup(self) -> None:
    with open(self.children_filename, 'wt', encoding='utf-8') as file:
      file.write('[]\n')
    self.children_memento = []
  


  def has_backed_up_children(self) -> bool:
    if not os.path.isfile(self.children_filename):
      return False
    with open(self.children_filename, 'rt', encoding='utf-8') as file:
      memento = json.load(file)
    if len(memento) == 0:
      return False
    self.children_memento = memento
    return True



  def restore_children_from_backup(self) -> List[Individual]:
    return [ Individual(**params) for params in self.children_memento ]



class Algorithm(ABC):
  """
  La representación genérica de algún algoritmo evolutivo. 
  """
  
  def __init__(self,
               workspaceName: str,
               referencePdbFilename: str,
               predictor: Predictor) -> None:
    """
    Constructor.
    - `referencePdbFilename`: el nombre del archivo PDB cuya estructura será 
      utilizada como referencia para calcular la función objetivo (o la función 
      de aptitud).
    - `predictor`: el algoritmo de predicción de la estructura de una proteína 
      a utilizar; este se utiliza para calcular la función objetivo (o la 
      función de aptitud).
    """
    super().__init__()
    self._predictor = predictor
    reference = PDB.load_structure_from_pdb_file(referencePdbFilename)
    self._sequence_length = Chain.count_residues_from_chain(reference)
    self._reference_backbone = Chain.filter_backbone_atoms_from_chain(reference)
    self.workspace = _Workspace(workspaceName, referencePdbFilename)



  @classmethod
  @abstractmethod
  def get_type(cls) -> str:
    pass



  @abstractmethod
  def run(self, iterationId: int, population: List[Individual]) -> None:
    pass



  def get_params_memento(self) -> dict:
    return {
      'workspaceName': self.workspace.name,
      'referencePdbFilename': self.workspace.reference_filename,
      'predictor': self._predictor.get_name()
    }



  def create_memento(self) -> dict:
    return {
      'algorithmType': self.get_type(),
      'algorithmParams': self.get_params_memento(),
      '__savedPopulations': self.workspace.saved_populations
    }
  


  def append_memento_and_cache(self, 
                               populations: Optional[List[str]] = None) -> None:
    if populations == None: populations = []
    self.workspace.saved_populations = populations
    self.workspace.memento = self.create_memento()



  def _get_pdb_filename(self, individual: Individual) -> str:
    return os.path.join(self.workspace.pdbs_folder, f'prot_{individual.id}.pdb')
