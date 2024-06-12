from .Metric import Metric
from ..Workspace import Workspace
import numpy as np
import numpy.typing as npt
import os





class ESM2Descriptors(Metric):

  _model = None
  _tokenizer = None



  @classmethod
  def _class_name(cls) -> str:
    return 'Metrics.ESM2Descriptors'
  


  def column_name(self) -> str:
    return ''
  


  def compute_value(self, **kwargs) -> None:
    sequence, sequence_id = kwargs['sequence'], kwargs['sequence_id']
    other_metrics = kwargs['otherMetrics']
    if self._model == None or self._tokenizer == None:
      from transformers import AutoTokenizer, EsmModel
      self._tokenizer = AutoTokenizer.from_pretrained('facebook/esm2_t36_3B_UR50D')
      self._model = EsmModel.from_pretrained('facebook/esm2_t36_3B_UR50D')
      # TODO usar este modelo en la GPU
      self._model.eval()
    workspace = Workspace.instance()
    csv_path = f'{workspace.esm2_dir}/{sequence_id}.csv'
    self.compute_descriptor_vector(sequence, csv_path)
    other_metrics['esm2_descriptor'] = csv_path
  


  def load_vector_from_csv(self, csvPath: str) -> npt.NDArray[np.float64]:
    for line in open(csvPath, 'rt'):
      values = line.split(',')
    vector = np.array([ float(x.strip()) for x in values ])
    return vector
  


  def compute_descriptor_vector(self, 
                                sequence: str, 
                                csvPath: str
                                ) -> npt.NDArray[np.float64]:
    workspace = Workspace.instance()
    os.makedirs(workspace.esm2_dir, exist_ok=True)
    import torch
    inputs = self._tokenizer(sequence, 
                             return_tensors='pt', 
                             padding=True, 
                             truncation=True)
    with torch.no_grad():
      outputs = self._model(**inputs)
    vector = outputs.last_hidden_state.numpy()
    vector = np.concatenate(vector, axis=1)
    vector = np.concatenate(vector, axis=0)
    with open(csvPath, 'wt', encoding='utf-8') as csv_file:
      values = vector.tolist()
      csv_file.write(','.join(values) + '\n')
    return vector
  