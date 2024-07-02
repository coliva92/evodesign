from .Metric import Metric
from .ESM2Descriptors import ESM2Descriptors
from ..Workspace import Workspace
import numpy as np





class ESM2DescriptorsRms(Metric):

  def column_name(self) -> str:
    return 'esm2_rms'
  


  def _params(self) -> dict:
    return { 'descriptor_calc': self._descriptor_calc.settings() }
  


  def __init__(self, descriptor_calc: ESM2Descriptors) -> None:
    super().__init__()
    self._ref_vector = None
    self._descriptor_calc = descriptor_calc 
  


  def compute_value(self, **kwargs) -> float:
    if self._ref_vector is None:
      ref_seq = kwargs['refSequence']
      workspace = Workspace.instance()
      ref_csv_path = f'{workspace.esm2_dir}/reference.csv'
      self._ref_vector = \
        self._descriptor_calc.compute_descriptor_vector(ref_seq, ref_csv_path)
    other_metrics = kwargs['otherMetrics']
    self._descriptor_calc(**kwargs)
    csv_path = other_metrics['esm2_descriptor']
    vector = self._descriptor_calc.load_vector_from_csv(csv_path)
    return np.sqrt(np.mean((self._ref_vector - vector) ** 2))
