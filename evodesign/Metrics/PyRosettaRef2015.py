from .Metric import Metric
from typing import Optional, List, Dict
from ..Context import Context
from Bio.PDB.Atom import Atom
import pandas as pd





class PyRosettaRef2015(Metric):

  score_fn = None
  


  def __init__(self, column: Optional[str] = None) -> None:
    super().__init__(column)
    if self.score_fn is not None: return
    import pyrosetta
    pyrosetta.init()
    self.score_fn = pyrosetta.get_score_function(True)
  


  def _compute_values(self, 
                      backbone: List[Atom],
                      data: pd.Series,
                      context: Context
                      ) -> pd.Series:
    pdb_path = f'{context.workspace.pdbs_dir}/prot_{data['sequence_id']}.pdb'
    import pyrosetta
    pose = pyrosetta.pose_from_pdb(pdb_path)
    score = self.score_fn(pose)
    data[self.column_name()] = score
    return data
  