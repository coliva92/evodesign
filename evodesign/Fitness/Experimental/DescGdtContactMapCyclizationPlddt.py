from typing import List
from .GdtContactMapCyclizationPlddt import GdtContactMapCyclizationPlddt
from ...Metrics.Experimental.iLearn.BLOSUM62 import BLOSUM62
from ...Metrics.Experimental.iLearn.CKSAAGP import CKSAAGP
from ...Metrics.Experimental.iLearn.CTDC import CTDC
from ...Metrics.Experimental.iLearn.CTDD import CTDD
from ...Metrics.Experimental.iLearn.CTDT import CTDT
from ...Metrics.Experimental.iLearn.EGAAC import EGAAC
from ...Metrics.Experimental.iLearn.GAAC import GAAC
from ...Metrics.Experimental.iLearn.GDPC import GDPC
from ...Metrics.Experimental.iLearn.GTPC import GTPC
import evodesign.Fitness.Experimental.Utils as Utils
import numpy as np
import pandas as pd





class DescGdtContactMapCyclizationPlddt(GdtContactMapCyclizationPlddt):

  _ref_desc = None



  @classmethod
  def _class_name(cls) -> str:
    return 'Fitness.Experimental.DescriptorsGdtContactMapCyclizationPlddt'
  


  @classmethod
  def column_name(cls) -> str:
    return 'fitness_desc_gdt_cm_cyc_plddt'
  


  def _params(self) -> dict:
    params = super()._params()
    params['ilearnDir'] = self._ilearn_dir
    params['referenceSequence'] = self._ref_sequence
    params['descWeight'] = self._desc_weight
    return params
  


  def __init__(self, 
               ilearnDir: str,
               referenceSequence: str,
               upperBound: float = 1.0, 
               cutoffs: List[float] = [ 1.0, 2.0, 4.0, 8.0 ], 
               gdtWeight: float = 1.0, 
               cmWeight: float = 1.0, 
               cycWeight: float = 1.0, 
               plddtWeight: float = 1.0,
               descWeight: float = 1.0
               ) -> None:
    super().__init__(upperBound, 
                     cutoffs, 
                     gdtWeight, 
                     cmWeight, 
                     cycWeight, 
                     plddtWeight)
    descriptors = [ 
      BLOSUM62(ilearnDir), CKSAAGP(ilearnDir), CTDC(ilearnDir), CTDD(ilearnDir),
      CTDT(ilearnDir), EGAAC(ilearnDir), GAAC(ilearnDir), GDPC(ilearnDir), 
      GTPC(ilearnDir)
    ]
    self._metrics = descriptors + self._metrics
    self._ilearn_dir = ilearnDir
    self._ref_sequence = referenceSequence
    self._desc_weight = descWeight
    self._weights = np.array([ descWeight, gdtWeight, cmWeight, cycWeight, plddtWeight ])
  


  def compute_fitness(self, **kwargs) -> float:
    if self._ref_desc == None: 
      self._compute_reference_descriptors()
    cyc = Utils.normalize(abs(Utils.cyclization_z_score(kwargs['cyclization'])))
    plddt = Utils.normalize_plddt(kwargs['plddt'])
    cm_rms = Utils.normalize(kwargs['contact_maps_rms'])
    temp = [
      pd.read_csv(kwargs['blosum62'], header=None).iloc[0][1:].to_numpy(),
      pd.read_csv(kwargs['CKSAAGP'], header=None).iloc[0][1:].to_numpy(),
      pd.read_csv(kwargs['CTDC'], header=None).iloc[0][1:].to_numpy(),
      pd.read_csv(kwargs['CTDD'], header=None).iloc[0][1:].to_numpy(),
      pd.read_csv(kwargs['CTDC'], header=None).iloc[0][1:].to_numpy(),
      pd.read_csv(kwargs['EGAAC'], header=None).iloc[0][1:].to_numpy(),
      pd.read_csv(kwargs['GAAC'], header=None).iloc[0][1:].to_numpy(),
      pd.read_csv(kwargs['GDPC'], header=None).iloc[0][1:].to_numpy(),
      pd.read_csv(kwargs['GTPC'], header=None).iloc[0][1:].to_numpy()
    ]
    temp[0] = Utils.normalize_numpy(temp[0], -4, 9)
    temp[3] = Utils.normalize_numpy(temp[3], 0, 100)
    model_desc = np.concatenate(temp)
    desc_rms = np.sqrt(np.mean((model_desc - self._ref_desc)**2))
    terms = np.array([ np.mean(desc_rms), kwargs['gdt'], cm_rms, cyc, plddt ])
    fitness = np.average(terms, weights=self._weights)
    return fitness
  


  def _compute_reference_descriptors(self):
    if self._ref_desc: return
    csv_paths = [
      metric(sequence=self._ref_sequence, sequenceId='ref_prot')
      for metric in self._metrics[0:8]
    ]
    temp = [
      pd.read_csv(csv_path, header=None).iloc[0][1:].to_numpy()
      for csv_path in csv_paths
    ]
    temp[0] = Utils.normalize_numpy(temp[0], -4, 11)
    temp[3] = Utils.normalize_numpy(temp[3], 0, 100)
    self._ref_desc = np.concatenate(temp)
