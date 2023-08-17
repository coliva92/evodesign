from abc import ABC, abstractmethod
from typing import List, Optional
from .Workspace import Workspace
from .Population import Population
import evodesign.Chain as Chain
from .Prediction import Predictor
from .Fitness import FitnessFunction
import random
import time





class Algorithm(ABC):

  @classmethod
  @abstractmethod
  def get_name(cls) -> str:
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
    reference = Chain.load_structure_from_pdb(targetPdbFilename)
    self._sequence_length = Chain.count_chain_residues(reference)
    self._reference_backbone = Chain.filter_backbone_atoms_in_chain(reference)
    self.best_solution = None
    self.workspace = Workspace(workspaceRoot,
                               targetPdbFilename)
    self.workspace.algorithm_settings = self.as_json()



  def as_json(self) -> dict:
    return {
      'algorithmType': self.get_name(),
      'algorithmParams': self._get_params_json()
    }
  


  def _get_params_json(self) -> dict:
    return {
      'workspaceRoot': self.workspace.root_folder,
      'targetPdbFilename': self.workspace.reference_filename,
      'predictor': self._predictor.get_name(),
      'fitnessFunction': self._fitness_fn.get_name()
    }
  


  @abstractmethod
  def __call__(self, population: Optional[Population] = None) -> None:
    raise NotImplementedError
