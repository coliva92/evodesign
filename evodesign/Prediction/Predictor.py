from abc import ABC, abstractmethod
from typing import List, Tuple
from Bio.PDB.Atom import Atom
import evodesign.Chain as Chain
from ..SettingsRetrievable import SettingsRetrievable
import numpy as np
import os





class Predictor(SettingsRetrievable, ABC):

  @abstractmethod
  def predict_structure(self, 
                        sequence: str, 
                        pdbPath: str
                        ) -> None:
    raise NotImplementedError
  


  @abstractmethod
  def predict_structure(self, 
                        sequence: str, 
                        pdbPath: str) -> None:
    raise NotImplementedError
  


  def __call__(self,
               sequence: str, 
               pdbPath: str
               ) -> Tuple[List[Atom], float]:
    # PDBParser cannot produce an instance of `Structure` directly from
    # a raw PDB string, it can only do it by reading from a PDB file.
    # Thus, the predicted structure must be stored first, before 
    # loading the `Structure` instance.
    if not os.path.isfile(pdbPath):
      os.makedirs(os.path.dirname(os.path.abspath(pdbPath)), exist_ok=True)
      self.predict_structure(sequence, pdbPath)
    structure = Chain.load_structure(pdbPath)
    bfactors = np.array([
      atom.get_bfactor() / 100.0 
        if atom.get_bfactor() > 1.0 
        else atom.get_bfactor()
      for atom in structure.get_atoms()
    ])
    return Chain.backbone_atoms(structure), bfactors.mean()
