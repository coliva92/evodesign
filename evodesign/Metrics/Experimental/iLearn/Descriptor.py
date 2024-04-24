from ...Metric import Metric
from abc import ABC
import subprocess
from ....Workspace import Workspace
import os





class Descriptor(ABC, Metric):

  def __init__(self, 
               ilearnDir: str,
               method: str
               ) -> None:
    self._ilearn_dir = ilearnDir
    self._method = method
  


  def __call__(self, **kwargs) -> str:
    sequence = kwargs['sequence']
    sequence_id = kwargs['sequence_id']
    workspace = Workspace.instance()
    fasta_path = f'{workspace.root_dir}/.temp.fasta'
    with open(fasta_path, 'wt', encoding='utf-8') as fasta_file:
      fasta_file.write(f'>{sequence_id}|null|null|null\n{sequence}\n')
    os.makedirs(workspace.ilearn_dir, exist_ok=True)
    csv_path = f'{workspace.ilearn_dir}/{sequence_id}_{self._method}.csv'
    for line in self._run_descriptor_method(self._method, fasta_path, csv_path):
      print(line, end='')
    os.remove(fasta_path)
    return csv_path
  


  def _run_descriptor_method(self,
                             method: str,
                             inputFastaPath: str,
                             outputCsvPath: str
                             ):
    cmd = [
      'python3',
      f'{self._ilearn_dir}/iLearn-protein-basic.py',
      '--method', method,
      '--format', 'csv',
      '--file', inputFastaPath,
      '--out', outputCsvPath
    ]
    popen = subprocess.Popen(cmd, 
                             stdout=subprocess.PIPE,
                             universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ''):
      yield stdout_line
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
      raise subprocess.CalledProcessError(return_code, cmd)
  