from .Metric import Metric
from ..Workspace import Workspace
from typing import Optional





class PyRosettaRef2015(Metric):

  score_fn = None
  


  def __init__(self, column: Optional[str] = None) -> None:
    super().__init__(column)
    if self.score_fn is not None: return
    import pyrosetta
    pyrosetta.init()
    self.score_fn = pyrosetta.get_score_function(True)
  


  def compute_value(self, **kwargs) -> float:
    workspace = Workspace.instance()
    pdb_path = f'{workspace.pdbs_dir}/prot_{kwargs["sequence"]}.pdb'
    import pyrosetta
    pose = pyrosetta.pose_from_pdb(pdb_path)
    return self.score_fn(pose)
  