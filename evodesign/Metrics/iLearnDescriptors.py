from .Metric import Metric
import subprocess
from ..Workspace import Workspace
from typing import Dict, List, Optional
import numpy.typing as npt
import numpy as np
import os





class iLearnDescriptors(Metric):
  
  def _params(self) -> dict:
    params = super()._params()
    params['ilearn_dir'] = self._ilearn_dir
    params['methods'] = self._methods
    return params



  def __init__(self, 
               ilearn_dir: str,
               methods: List[str] = [
                 'GAAC', 'EGAAC', 'CKSAAGP', 'GDPC', 'GTPC', 'BLOSUM62', 'CTDC',
                 'CTDT', 'CTDD'
               ],
               column: Optional[str] = None
               ) -> None:
    super().__init__(column)
    self._ilearn_dir = ilearn_dir
    self._methods = methods
  


  def compute_value(self, **kwargs) -> None:
    other_metrics = kwargs['otherMetrics']
    sequence = kwargs['sequence']
    sequence_id = kwargs['sequenceId']
    workspace = Workspace.instance()
    fasta_path = f'{workspace.root_dir}/.temp.fasta'
    with open(fasta_path, 'wt', encoding='utf-8') as fasta_file:
      fasta_file.write(f'>{sequence_id}|null|null|null\n{sequence}\n')
    csv_dir = self.vectors_csv_dir(sequence_id)
    vector_paths = self.compute_descriptors(fasta_path, csv_dir)
    for method, filepath in vector_paths.items():
      other_metrics[method] = filepath
    os.remove(fasta_path)
  


  def vectors_csv_dir(self, sequence_id: str) -> str:
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
    vector = np.array([ float(x.strip()) for x in values if x.find('null') == -1 ])
    return vector
  


  def compute_descriptors(self, 
                          fastaPath: str,
                          csvDir: str
                          ) -> Dict[str, str]:
    os.makedirs(csvDir, exist_ok=True)
    output = {}
    for method in self._methods:
      csv_path = f'{csvDir}/{method}.csv'
      output[method] = csv_path
      for line in self._descriptor_method_cmd(method, fastaPath, csv_path):
        print(line, end='')
    return output



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
  