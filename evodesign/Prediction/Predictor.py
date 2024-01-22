from abc import ABC, abstractmethod, abstractclassmethod
from typing import List
from Bio.PDB.Atom import Atom
from Bio.PDB import PDBParser
from evodesign import Chain
import os
import time





class Predictor(ABC):

  _pdb_parser = None



  @abstractclassmethod
  def name(cls) -> str:
    raise NotImplementedError


  
  @abstractmethod
  def predict_structure(self, 
                        sequence: str, 
                        pdbFilename: str
                        ) -> None:
    raise NotImplementedError
  


  def __call__(self,
               sequence: str, 
               pdbFilename: str
               ) -> List[Atom]:
    # PDBParser cannot produce an instance of `Structure` directly from
    # a raw PDB string, it can only do it by reading from a PDB file.
    # Thus, the predicted structure must be stored first, before 
    # loading the `Structure` instance.
    if not os.path.isfile(pdbFilename):
      os.makedirs(os.path.dirname(os.path.abspath(pdbFilename)), exist_ok=True)
      self.predict_structure(sequence, pdbFilename)
    if not self._pdb_parser:
      self._pdb_parser = PDBParser()
    structure = self._pdb_parser.get_structure(str(time.time_ns()), pdbFilename)
    return Chain.backbone_atoms(structure)
