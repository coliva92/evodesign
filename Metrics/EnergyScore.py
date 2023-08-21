from .Metric import Metric
from typing import List
from Bio.PDB.Atom import Atom
import subprocess
import os





class EnergyScore(Metric):
  
  def __init__(self) -> None:
    super().__init__()
    import pyrosetta
    pyrosetta.init()
    from pyrosetta.teaching import get_score_function
    self._score_fn = get_score_function(True)
  


  def initialize(self, 
                 referencePdbFilename: str,
                 workspaceRoot: str,
                 pdbsFolder: str
                 ) -> None:
    self._pdbs_folder = pdbsFolder
    self._sequence_filename = os.path.join(workspaceRoot, 
                                           '~sequence.tmp')
    self._scwrl_command = [
      os.path.join('./scwrl', 'Scwrl4'),
      '-i', referencePdbFilename,
      '-o', None,
      '-s', self._sequence_filename,
      '-h', '-t'
    ]

  

  def __call__(self, 
               modelBackbone: List[Atom], 
               referenceBackbone: List[Atom],
               sequence: str
               ) -> float:
    with open(self._sequence_filename, 'wt', encoding='utf-8') as seq_file:
      seq_file.write(f'{sequence}\n')
    side_chain_packing_filename = os.path.join(self._pdbs_folder, 
                                               f'scwrl_{sequence}.pdb')
    self._scwrl_command[4] = side_chain_packing_filename
    subprocess.run(self._scwrl_command,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL)
    import pyrosetta
    pose = pyrosetta.pose_from_pdb(side_chain_packing_filename)
    return self._score_fn(pose)
    