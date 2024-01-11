from . import Metric
from typing import List, Optional
from Bio.PDB.Atom import Atom
from Bio.PDB import Superimposer





class Rmsd(Metric):
  
  def __call__(self, 
               modelBackbone: List[Atom], 
               referenceBackbone: List[Atom],
               _: Optional[str] = None
               ) -> float:
    superimposer = Superimposer()
    superimposer.set_atoms(referenceBackbone, modelBackbone)
    superimposer.apply(modelBackbone)
    return superimposer.rms
  