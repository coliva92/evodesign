from abc import ABC, abstractmethod, abstractclassmethod
from typing import List
from Bio.PDB.Atom import Atom
from ..Chain import Chain
from ..AsSetttings import AsSettings
import os





class Predictor(AsSettings, ABC):

  @abstractmethod
  def predict_structure(self, 
                        sequence: str, 
                        pdbPath: str
                        ) -> None:
    raise NotImplementedError
  


  def __call__(self,
               sequence: str, 
               pdbPath: str
               ) -> List[Atom]:
    # PDBParser cannot produce an instance of `Structure` directly from
    # a raw PDB string, it can only do it by reading from a PDB file.
    # Thus, the predicted structure must be stored first, before 
    # loading the `Structure` instance.
    if not os.path.isfile(pdbPath):
      os.makedirs(os.path.dirname(os.path.abspath(pdbPath)), exist_ok=True)
      self.predict_structure(sequence, pdbPath)
    structure = Chain.load_structure(pdbPath)
    return Chain.backbone_atoms(structure)
