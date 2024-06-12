from ..Metric import Metric
from ..iLearnDescriptors import iLearnDescriptors
from ..ESM2Descriptors import ESM2Descriptors
from ...Workspace import Workspace
from typing import List
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
  _ref_vectors = None



  @classmethod
  def _class_name(cls) -> str:
    return 'Metrics.Yan24.Descriptors'
  


  def column_name(self) -> str:
    return 'yan24_descriptors'
  


  def _params(self) -> dict:
    return {
      'weight': self._weight,
      'ilearnDir': self._ilearn_dir,
      'ilearnMethods': self._ilearn_methods
    }
  


  def __init__(self,
               weight: float,
               ilearnDir: str,
               ilearnMethods: List[str]
               ) -> None:
    super().__init__()
    self._ilearn_model = iLearnDescriptors(ilearnDir, ilearnMethods)
    self._esm2_model = ESM2Descriptors()
    self._ilearn_dir = ilearnDir
    self._ilearn_methods = ilearnMethods
    self._weight = weight
  


  def compute_value(self, **kwargs) -> float:
    if self._ref_vectors == None:
      workspace = Workspace.instance()
      self._ilearn_model.compute_descriptors(workspace.target_fasta_path, 
                                             workspace.ilearn_dir)
      self._ref_vectors = self._ilearn_model.load_vectors(workspace.ilearn_dir)
      csv_path = f'{workspace.esm2_dir}/reference.csv'
      self._ref_vectors['ESM2'] = \
        self._esm2_model.compute_descriptor_vector(kwargs['refSequence'], 
                                                   csv_path)
    self._ilearn_model(**kwargs)
    vectors_dir = self._ilearn_model.vectors_csv_path(kwargs['sequence_id'])
    model_vectors = self._ilearn_model.load_vectors(vectors_dir)
    csv_path = f'{workspace.esm2_dir}/{kwargs["sequence_id"]}.csv'
    model_vectors['EMS2'] = \
      self._esm2_model.compute_descriptor_vector(kwargs['sequence'], csv_path)
    divergence_scores = np.array([
      self._symmetric_kullback_leibler(self._ref_vectors, model_vectors, method)
      for method in self._ilearn_methods
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
    if methodType == 'non_probabilistic':
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
    divergence = sp.entropy(a, b)
    return divergence
