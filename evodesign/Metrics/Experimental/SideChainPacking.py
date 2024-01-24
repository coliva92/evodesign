from ..Metric import Metric
from typing import List
from Bio.PDB.Atom import Atom
from ...Workspace import Workspace
import subprocess
import os





class SideChainPacking(Metric):

  def __init__(self,
               scwrlExecutablePath: str = './scwrl4/Scwrl4'
               ) -> None:
    super().__init__()
    workspace = Workspace.instance()
    self._scwrl_outputs_folder = f'{workspace.root_dir}/scwrl_outputs'
    self._scwrl_sequence_filename = f'{workspace.root_dir}/.scwrl4_in'
    self._scwrl_input_pdb_filename = f'{workspace.target_pdb_path}'
    self._scwrl_executable = scwrlExecutablePath



  def __call__(self, 
               sequence: str,
               model: List[Atom], 
               reference: List[Atom]
               ) -> float:
    os.makedirs(self._scwrl_outputs_folder, exist_ok=True)
    with open(self._scwrl_sequence_filename, 'wt', encoding='utf-8') as seq_file:
      seq_file.write(f'{sequence}\n')
    output_pdb_filename = f'{self._scwrl_outputs_folder}/prot_{sequence}.pdb'
    output_filename = f'{self._scwrl_outputs_folder}/scwrl4_output_{sequence}.txt'
    with open(output_filename, 'wt', encoding='utf-8') as output_file:
      args = [
        self._scwrl_executable,
        '-i', self._scwrl_input_pdb_filename,
        '-s', self._scwrl_sequence_filename,
        '-o', output_pdb_filename,
        '-h', # do not write H atoms in the resulting PDB file
        '-t'  # do not write the OXT atom in the resulting PDB file
      ]
      subprocess.run(args,
                     stdout=output_file,
                     stderr=subprocess.DEVNULL)
    for line in open(output_filename, 'r'):
      if line.find('Total minimal energy of the graph') != -1:
        return float(line[36:].strip())
