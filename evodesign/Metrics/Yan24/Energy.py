from ..Metric import Metric
from ..PyRosettaRef2015 import PyRosettaRef2015
from .Geometric import Geometric
from ...Workspace import Workspace
from typing import List
import numpy as np





class Energy(Metric):

  _reference_energy = None



  @classmethod
  def _class_name(cls) -> str:
    return 'Metrics.Yan24.Energy'



  def column_name(self) -> str:
    return 'yan24_energy'
  


  def _params(self) -> dict:
    return { 
      'weight': self._weight,
      'cutoffs': self._cutoffs
    }
  


  def __init__(self, 
               weight: float,
               cutoffs: List[float] = [ 1.0, 2.0, 4.0, 8.0 ]
               ) -> None:
    super().__init__()
    self._weight = weight
    self._cutoffs = cutoffs
    self._energy_calc = PyRosettaRef2015()
    self._geometric_calc = Geometric(cutoffs)
  


  def compute_value(self, **kwargs) -> float:
    if self._reference_energy == None:
      workspace = Workspace.instance()
      import pyrosetta
      pose = pyrosetta.pose_from_pdb(workspace.target_pdb_path)
      self._reference_energy = self._energy_calc.score_fn(pose)
    other_metrics = kwargs['otherMetrics']
    e = other_metrics[self._energy_calc.column_name()] = \
      self._energy_calc(**kwargs)
    f1 = other_metrics[self._geometric_calc.column_name()] = \
      self._geometric_calc(**kwargs)
    return f1 / 4 + 1 / np.exp(abs(self._reference_energy - e) / self._weight)
  