from .Metric import Metric
from .ESM2Descriptors import ESM2Descriptors
from ..Workspace import Workspace
import numpy as np





class ESM2DescriptorsRms(Metric):

  def column_name(self) -> str:
    return 'esm2_rms'
  


  def _params(self) -> dict:
    return self._esm2_desc._params()
  


  def __init__(self, use_gpu: bool = True) -> None:
    super().__init__()
    self._ref_vectors = None
    self._esm2_desc = ESM2Descriptors(use_gpu)
  


  def compute_value(self, **kwargs) -> float | None:
    if not self._ref_vectors:
      ref_seq = kwargs['refSequence']
      workspace = Workspace.instance()
      ref_csv_path = f'{workspace.esm2_dir}/reference.csv'
      self._ref_vectors = \
        self._esm2_desc.compute_descriptor_vector(ref_seq, ref_csv_path)
    other_metrics = kwargs['otherMetrics']
    self._esm2_desc(**kwargs)
    csv_path = other_metrics['esm2_descriptor']
    vectors = self._esm2_desc.load_vector_from_csv(csv_path)
    temp = [
      np.sqrt(np.mean((u - v) ** 2))
      for u, v in zip(self._ref_vectors, vectors)
    ]
    return np.mean(temp)
