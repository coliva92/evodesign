from ..Metric import Metric
from ..PyRosettaRef2015 import PyRosettaRef2015
from .Geometry import Geometry
from ...Workspace import Workspace
from typing import List, Optional
import numpy as np





class Energy(Metric):

  _reference_energy = None
  


  def _params(self) -> dict:
    params = super()._params()
    params['weight'] = self._weight
    params['geometry_metric'] = self._geometry_metric.settings()
    return params
  


  def __init__(self, 
               weight: float,
               geometry_metric: Optional[Geometry] = None,
               column: Optional[str] = None
               ) -> None:
    super().__init__(column)
    self._weight = weight
    self._geometry_metric = geometry_metric
    if self._geometry_metric is None:
      self._geometry_metric = Geometry()
    self._energy_calc = PyRosettaRef2015()
  


  def compute_value(self, **kwargs) -> float:
    if self._reference_energy is None:
      workspace = Workspace.instance()
      import pyrosetta
      pose = pyrosetta.pose_from_pdb(workspace.target_pdb_path)
      self._reference_energy = self._energy_calc.score_fn(pose)
    other_metrics = kwargs['otherMetrics']
    e = other_metrics[self._energy_calc.column_name()] = \
      self._energy_calc(**kwargs)
    f1 = other_metrics[self._geometry_metric.column_name()] = \
      self._geometry_metric(**kwargs)
    return f1 / 4 + 1 / np.exp(abs(self._reference_energy - e) / self._weight)
  