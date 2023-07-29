from abc import ABC, abstractmethod
from typing import List, Optional
from .Workspace import Workspace
from .Population import Individual
import evodesign.Chain as Chain
from .Prediction import Predictor
import os





class Algorithm(ABC):
  """La representación genérica de algún algoritmo evolutivo. 
  """
  
  @classmethod
  @abstractmethod
  def get_name(cls) -> str:
    """Retorna el nombre del algoritmo.
    """
    raise NotImplementedError
  


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
    reference = Chain.load_structure_from_pdb(targetPdbFilename)
    self._sequence_length = Chain.count_chain_residues(reference)
    self._reference_backbone = Chain.filter_backbone_atoms_in_chain(reference)
    self.workspace = Workspace(workspaceName, targetPdbFilename)
    self.best_solution = None



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
      'algorithmType': self.get_name(),
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
