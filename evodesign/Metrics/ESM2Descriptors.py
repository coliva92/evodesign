from .Metric import Metric
from ..Workspace import Workspace
import numpy as np
import numpy.typing as npt
import os





class ESM2Descriptors(Metric):

  model = None
  batch_converter = None
  


  def column_name(self) -> str:
    return ''
  


  def _params(self) -> dict:
    return { 'useGpu': self._use_gpu }
  


  def __init__(self, useGpu: bool = True) -> None:
    self._use_gpu = useGpu
  


  def compute_value(self, **kwargs) -> None:
    sequence, sequence_id = kwargs['sequence'], kwargs['sequence_id']
    other_metrics = kwargs['otherMetrics']
    if self.model == None or self.batch_converter == None:
      import esm
      import torch
      self.model, alphabet = esm.pretrained.esm2_t36_3B_UR50D()
      self.batch_converter = alphabet.get_batch_converter()
      self.model.eval()
      if self._use_gpu and torch.cuda.is_available(): 
        self.model = self.model.cuda()
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
                                csvPath: str,
                                sequence_id: str = 'prot_0000_0000'
                                ) -> npt.NDArray[np.float64]:
    import torch
    workspace = Workspace.instance()
    os.makedirs(workspace.esm2_dir, exist_ok=True)
    data = [( sequence_id, sequence )]
    seq_ids, seqs, tokens = self.batch_converter(data)
    if self._use_gpu and torch.cuda.is_available():
      tokens = tokens.to(device='cuda', non_blocking=True)
    with torch.no_grad():
      result = self.model(tokens, 
                          repr_layers=[self.model.num_layers], 
                          return_contacts=False)
    
    # `result['representations']` contains the weights of each layer in the 
    # NN; we only want the weights of the last layer
    last_layer = result['representations'][self.model.num_layers]

    # the last layer contains a certain number of weights per token; a token is
    # an integer representation of each AA in the input sequence, however, 
    # additional tokens are appended at the beginning and at the end of said 
    # sequence; we only want to retrieve the weights corresponding to the AA
    # in the sequence
    matrix = last_layer[0][1 : len(seqs[0]) + 1].cpu().numpy()

    # free GPU memory
    del tokens
    del result
    torch.cuda.empty_cache()

    # collapse all weights into a single vector
    vector = np.concatenate(matrix, axis=0)

    # save the descriptors to a file
    with open(csvPath, 'wt', encoding='utf-8') as csv_file:
      values = vector.tolist()
      csv_file.write(','.join(values) + '\n')
    return vector
  