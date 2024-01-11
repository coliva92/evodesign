from abc import ABC, abstractmethod, abstractclassmethod
from typing import Optional
from Workspace import Workspace
from Population import Population
import Chain
import FileIO
from Prediction import Predictor
from Fitness import FitnessFunction





class Algorithm(ABC):

  @abstractclassmethod
  def name(cls) -> str:
    raise NotImplementedError
  


  def __init__(self,
               workspaceRoot: str,
               targetPdbFilename: str,
               predictor: Predictor,
               fitnessFunction: FitnessFunction
               ) -> None:
    super().__init__()
    self._predictor = predictor
    self._fitness_fn = fitnessFunction
    reference = FileIO.load_structure_from_pdb(targetPdbFilename)
    self._sequence_length = Chain.count_chain_residues(reference)
    self._reference_backbone = Chain.filter_backbone_atoms_in_chain(reference)
    self.best_solution = None
    self.workspace = Workspace(workspaceRoot, targetPdbFilename)



  def settings(self) -> dict:
    return {
      'algorithmType': self.name(),
      'algorithmParams': self.params_as_dict()
    }
  


  def params_as_dict(self) -> dict:
    return {
      'workspaceRoot': self.workspace.root_folder,
      'targetPdbFilename': self.workspace.reference_filename,
      'predictor': self._predictor.name(),
      'fitnessFunction': self._fitness_fn.name(),
      'fitnessFunctionParams': self._fitness_fn.params_as_dict()
    }
  


  @abstractmethod
  def __call__(self, population: Optional[Population] = None) -> None:
    raise NotImplementedError
