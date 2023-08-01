from .Metric import Metric
from typing import List
from Bio.PDB.Atom import Atom
from Bio.PDB import Superimposer





class Rmsd(Metric):
  
  def __call__(self, 
               modelBackbone: List[Atom], 
               referenceBackbone: List[Atom]
               ) -> float:
    superimposer = Superimposer()
    superimposer.set_atoms(referenceBackbone, modelBackbone)
    superimposer.apply(modelBackbone)
    return superimposer.rms
  