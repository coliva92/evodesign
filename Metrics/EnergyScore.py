from .Metric import Metric
from typing import List
from Bio.PDB.Atom import Atom





class EnergyScore(Metric):
  
  def __init__(self) -> None:
    super().__init__()
    import pyrosetta
    pyrosetta.init()
    from pyrosetta.teaching import get_score_function
    self._score_fn = get_score_function(True)

  

  def __call__(self, 
               modelBackbone: List[Atom], 
               referenceBackbone: List[Atom],
               sequence: str
               ) -> float:
    filename = f'test/pdbs/prot_{sequence}.pdb'
    import pyrosetta
    pose = pyrosetta.pose_from_pdb(filename)
    return self._score_fn(pose)
    