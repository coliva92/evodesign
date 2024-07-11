from ..Metric import Metric
from ..iLearnDescriptors import iLearnDescriptors
from ..ESM2Descriptors import ESM2Descriptors
from ...Workspace import Workspace
from typing import List, Optional
import numpy as np
import numpy.typing as npt
import scipy.stats as sp





class Descriptors(Metric):

  _DESCRIPTOR_METHOD_TYPES = {
    'GAAC': 'probabilistic',
    'EGAAC': 'non_probabilistic', 
    'CKSAAGP': 'non_probabilistic',
    'GDPC': 'probabilistic',
    'GTPC': 'probabilistic',
    'BLOSUM62': 'non_probabilistic',
    'CTDC': 'non_probabilistic',
    'CTDT': 'non_probabilistic',
    'CTDD': 'non_probabilistic',
    'ESM2': 'non_probabilistic'
  }
  


  def _params(self) -> dict:
    params = self._params()
    params['weight'] = self._weight
    params['ilearm_metric'] = self._ilearn_metric.settings()
    params['esm2_metric'] = self._esm2_metric.settings()
    return params
  


  def __init__(self,
               weight: float,
               ilearn_metric: iLearnDescriptors,
               esm2_metric: Optional[ESM2Descriptors] = None,
               column: Optional[str] = None
               ) -> None:
    super().__init__(column)
    self._weight = weight
    self._ilearn_metric = ilearn_metric
    if esm2_metric is None:
      esm2_metric = ESM2Descriptors()
    self._esm2_metric = esm2_metric
    self._ref_vectors = None
  


  def compute_value(self, **kwargs) -> float:
    workspace = Workspace.instance()
    if self._ref_vectors is None:
      self._ilearn_metric.compute_descriptors(workspace.target_fasta_path, 
                                             workspace.ilearn_dir)
      self._ref_vectors = self._ilearn_metric.load_vectors(workspace.ilearn_dir)
      csv_path = f'{workspace.esm2_dir}/reference.csv'
      self._ref_vectors['ESM2'] = \
        self._esm2_metric.compute_descriptor_vector(kwargs['refSequence'], 
                                                   csv_path)
    self._ilearn_metric(**kwargs)
    vectors_dir = self._ilearn_metric.vectors_csv_dir(kwargs['sequenceId'])
    model_vectors = self._ilearn_metric.load_vectors(vectors_dir)
    csv_path = f'{workspace.esm2_dir}/{kwargs["sequenceId"]}.csv'
    model_vectors['EMS2'] = \
      self._esm2_metric.compute_descriptor_vector(kwargs['sequence'], csv_path)
    divergence_scores = np.array([
      self._symmetric_kullback_leibler(self._ref_vectors[method], 
                                       model_vectors[method], 
                                       self._DESCRIPTOR_METHOD_TYPES[method])
      for method in self._ilearn_metric._methods
    ])
    return 2 * self._weight / (1 + divergence_scores.mean())
  


  def _symmetric_kullback_leibler(self,
                                  reference: npt.NDArray[np.float64],
                                  model: npt.NDArray[np.float64],
                                  methodType: str
                                  ) -> float:
    a = self._kullback_leibler(reference, model, methodType)
    b = self._kullback_leibler(model, reference, methodType)
    return a + b



  def _kullback_leibler(self,
                        reference: npt.NDArray[np.float64],
                        model: npt.NDArray[np.float64],
                        methodType: str
                        ) -> float:
    a, b = reference, model
    if methodType == 'probabilistic':
      # we do this so that the probability is non-zero
      a += 0.01
      b += 0.01
    else:
      min_value = min([ reference.min(), model.min() ])
      max_value = max([ reference.max(), model.max() ])
      dist_range = (min_value, max_value)
      a, _ = np.histogram(reference, 
                          bins=20, 
                          range=dist_range, 
                          density=False)
      b, _ = np.histogram(model, 
                          bins=20, 
                          range=dist_range,
                          density=False)
      # numpy does not like mixing int with float
      if a.dtype != np.float64: a = a.astype(np.float64)
      if b.dtype != np.float64: b = b.astype(np.float64)
      # we do this so that the probability is non-zero
      a += 0.5
      b += 0.5
    divergence = sp.entropy(a, b)
    return divergence
