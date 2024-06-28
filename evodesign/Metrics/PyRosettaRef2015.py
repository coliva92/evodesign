from .Metric import Metric
from ..Workspace import Workspace





class PyRosettaRef2015(Metric):

  score_fn = None
  

  
  def column_name(self) -> str:
    return 'ref2015'
  


  def __init__(self) -> None:
    super().__init__()
    if self.score_fn != None: return
    import pyrosetta
    pyrosetta.init()
    self.score_fn = pyrosetta.get_score_function(True)
  


  def compute_value(self, **kwargs) -> float:
    workspace = Workspace.instance()
    pdb_path = f'{workspace.pdbs_dir}/prot_{kwargs["sequence"]}.pdb'
    import pyrosetta
    pose = pyrosetta.pose_from_pdb(pdb_path)
    return self.score_fn(pose)
  