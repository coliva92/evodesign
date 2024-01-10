from .Metric import Metric
from typing import List, Optional
from Bio.PDB.Atom import Atom
import subprocess
import os





class SideChainPackingEnergyScore(Metric):

  def __init__(self,
               workspaceRoot: str,
               targetPdbFilename: str,
               scwrlExecutablePath: str = './scwrl4/Scwrl4'
               ) -> None:
    super().__init__()
    self._scwrl_outputs_folder = os.path.join(workspaceRoot, 'scwrl_outputs')
    self._scwrl_sequence_filename = os.path.join(workspaceRoot, '.scwrl4.in')
    self._scwrl_input_pdb_filename = targetPdbFilename
    self._scwrl_executable = scwrlExecutablePath
    self._root_folder = workspaceRoot



  def __call__(self, 
               modelBackbone: List[Atom], 
               referenceBackbone: List[Atom],
               sequence: Optional[str] = None
               ) -> float:
    os.makedirs(self._scwrl_outputs_folder, exist_ok=True)
    with open(self._scwrl_sequence_filename, 'wt', encoding='utf-8') as seq_file:
      seq_file.write(f'{sequence}\n')
    output_pdb_filename = os.path.join(self._scwrl_outputs_folder, 
                                       f'prot_{sequence}.pdb')
    output_filename = os.path.join(self._scwrl_outputs_folder,
                                   f'scwrl4_output_{sequence}.txt')
    with open(output_filename, 'wt', encoding='utf-8') as output_file:
      args = [
        self._scwrl_executable,
        '-i', self._scwrl_input_pdb_filename,
        '-s', self._scwrl_sequence_filename,
        '-o', output_pdb_filename,
        '-h', # no se escriben los átomos de H en el PDB resultante
        '-t'  # no se escribe el átomo OXT en el PDB resultante
      ]
      subprocess.run(args,
                     stdout=output_file,
                     stderr=subprocess.DEVNULL)
    for line in open(output_filename, 'r'):
      if line.find('Total minimal energy of the graph') != -1:
        return float(line[36:].strip())
