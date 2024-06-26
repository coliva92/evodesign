from .Metric import Metric
import subprocess
from ..Workspace import Workspace
from typing import Dict, List
import numpy.typing as npt
import numpy as np
import os





class iLearnDescriptors(Metric):
  
  def column_name(self) -> str:
    return ''
  


  def _params(self) -> dict:
    params = super()._params()
    params['ilearnDir'] = self._ilearn_dir
    params['methods'] = self._methods
    return params



  def __init__(self, 
               ilearnDir: str,
               methods: List[str]
               ) -> None:
    self._ilearn_dir = ilearnDir
    self._methods = methods
  


  def compute_value(self, **kwargs) -> None:
    other_metrics = kwargs['otherMetrics']
    sequence = kwargs['sequence']
    sequence_id = kwargs['sequenceId']
    workspace = Workspace.instance()
    fasta_path = f'{workspace.root_dir}/.temp.fasta'
    with open(fasta_path, 'wt', encoding='utf-8') as fasta_file:
      fasta_file.write(f'>{sequence_id}|null|null|null\n{sequence}\n')
    csv_path = self.vectors_csv_path(sequence_id)
    vector_paths = self.compute_descriptors(fasta_path, csv_path)
    for method, filepath in vector_paths.items():
      other_metrics[method] = filepath
    os.remove(fasta_path)
  


  def vectors_csv_path(self, sequence_id: str) -> str:
    workspace = Workspace.instance()
    return f'{workspace.ilearn_dir}/{sequence_id}'
  

  
  def load_vectors(self, vectorsDir: str) -> Dict[str, npt.NDArray[np.float64]]:
    descriptors = {
      method: self.load_vector_from_csv(f'{vectorsDir}/{method}.csv')
      for method in self._methods
    }
    return descriptors
  


  def load_vector_from_csv(self, csvPath: str) -> npt.NDArray[np.float64]:
    for line in open(csvPath, 'rt'):
      values = line.split(',')
    vector = np.array([ float(x.strip()) for x in values ])
    return vector
  


  def compute_descriptors(self, 
                          fastaPath: str,
                          csvDir: str
                          ) -> Dict[str, str]:
    workspace = Workspace.instance()
    os.makedirs(workspace.ilearn_dir, exist_ok=True)
    output = {}
    for method in self._methods:
      csv_path = f'{csvDir}/{method}.csv'
      output[method] = csv_path
      for line in self._descriptor_method_cmd(method, fastaPath, csv_path):
        print(line, end='')
    return



  def _descriptor_method_cmd(self,
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
  