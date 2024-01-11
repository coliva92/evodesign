from abc import ABC, abstractmethod, abstractclassmethod
from typing import List
from Bio.PDB.Atom import Atom
from Bio.PDB import PDBParser
import Chain
import os
import time





class Predictor(ABC):

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
    if not os.path.isfile(pdbFilename):
      os.makedirs(os.path.dirname(os.path.abspath(pdbFilename)), exist_ok=True)
      # PDBParser no puede producir una instancia de la clase `Structure` 
      # directamente de datos crudos, solo puede hacerlo leyendo directamente 
      # un achivo PDB; entonces, la estructura predicha debe guardarse en un 
      # archivo PDB
      self.predict_structure(sequence, pdbFilename)
    parser = PDBParser()
    structure = parser.get_structure(time.time_ns(), pdbFilename)
    return Chain.filter_backbone_atoms_in_chain(structure)
