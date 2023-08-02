from abc import ABC, abstractmethod
from typing import List, Optional
from .Workspace import Workspace
from .Individual import Individual
import evodesign.Chain as Chain
from .Prediction import Predictor
from .Fitness import FitnessFunction





class Algorithm(ABC):
  
  @classmethod
  @abstractmethod
  def get_name(cls) -> str:
    raise NotImplementedError
  


  def __init__(self,
               workspaceName: str,
               targetPdbFilename: str,
               predictor: Predictor,
               fitnessFunction: FitnessFunction,
               populationFilenames: Optional[List[str]] = None
               ) -> None:
    super().__init__()
    if populationFilenames is None: populationFilenames = []
    self._predictor = predictor
    self._fitness_fn = fitnessFunction
    reference = Chain.load_structure_from_pdb(targetPdbFilename)
    self._sequence_length = Chain.count_chain_residues(reference)
    self._reference_backbone = Chain.filter_backbone_atoms_in_chain(reference)
    self.workspace = Workspace(workspaceName, 
                               targetPdbFilename, 
                               populationFilenames)
    self.workspace.memento = self.create_memento()
    self.best_solution = None



  @abstractmethod
  def run(self, 
          iterationId: int, 
          population: List[Individual]
          ) -> None:
    raise NotImplementedError



  def create_memento(self) -> dict:
    return {
      'algorithmType': self.get_name(),
      'algorithmParams': self._get_params_memento(),
      '__savedPopulations': self.workspace.population_filenames
    }
  


  def _get_params_memento(self) -> dict:
    return {
      'workspaceName': self.workspace.name,
      'targetPdbFilename': self.workspace.reference_filename,
      'predictor': self._predictor.get_name(),
      'fitnessFunction': self._fitness_fn.get_name()
    }
