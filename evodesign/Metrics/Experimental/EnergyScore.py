from ..Metric import Metric
from typing import List, Optional
from Bio.PDB.Atom import Atom
from ...Workspace import Workspace





class EnergyScore(Metric):

  score_fn = None



  @classmethod
  def column_name(cls) -> str:
    return 'Energy_Score'
  

  
  def __init__(self) -> None:
    super().__init__()
    if not self.score_fn:
      import pyrosetta
      pyrosetta.init()
      from pyrosetta.teaching import get_score_function
      self.score_fn = get_score_function(True)
  


  def __call__(self, 
               model: List[Atom], 
               reference: List[Atom],
               sequence: Optional[str] = None
               ) -> float:
    workspace = Workspace.instance()
    filename = f'{workspace.pdbs_dir}/prot_{sequence}.pdb'
    import pyrosetta
    pose = pyrosetta.pose_from_pdb(filename)
    return self.score_fn(pose)
    