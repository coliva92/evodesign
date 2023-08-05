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
               workspaceName: str,
               targetPdbFilename: str,
               predictor: Predictor,
               fitnessFunction: FitnessFunction,
               populationFilenames: Optional[List[str]] = None,
               seed: Optional[int] = None
               ) -> None:
    super().__init__()
    if populationFilenames is None: populationFilenames = []
    if seed is None: seed = time.time()
    random.seed(seed)
    self._seed = seed
    self._predictor = predictor
    self._fitness_fn = fitnessFunction
    reference = Chain.load_structure_from_pdb(targetPdbFilename)
    self._sequence_length = Chain.count_chain_residues(reference)
    self._reference_backbone = Chain.filter_backbone_atoms_in_chain(reference)
    self.workspace = Workspace(workspaceName, 
                               self.as_json(),
                               targetPdbFilename, 
                               populationFilenames)
    self.best_solution = None



  def as_json(self) -> dict:
    return {
      'algorithmType': self.get_name(),
      'algorithmParams': self._get_params_json(),
    }
  


  def _get_params_json(self) -> dict:
    return {
      'workspaceName': self.workspace.name,
      'seed': self._seed,
      'targetPdbFilename': self.workspace.reference_filename,
      'predictor': self._predictor.get_name(),
      'fitnessFunction': self._fitness_fn.get_name()
    }
  


  @abstractmethod
  def __call__(self, population: Optional[Population] = None) -> None:
    raise NotImplementedError
