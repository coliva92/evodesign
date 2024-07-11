from ..Metric import Metric
from ..PyRosettaRef2015 import PyRosettaRef2015
from .Geometry import Geometry
from ...Workspace import Workspace
from typing import Optional
import numpy as np





class Energy(Metric):
  
  def _params(self) -> dict:
    params = super()._params()
    params['weight'] = self._weight
    params['geometry_metric'] = self._geometry_metric.settings()
    params['energy_metric'] = self._energy_metric
    return params
  


  def __init__(self, 
               weight: float,
               geometry_metric: Optional[Geometry] = None,
               energy_metric: Optional[PyRosettaRef2015] = None,
               column: Optional[str] = None
               ) -> None:
    super().__init__(column)
    self._weight = weight
    if geometry_metric is None:
      geometry_metric = Geometry()
    self._geometry_metric = geometry_metric
    if energy_metric is None:
      energy_metric = PyRosettaRef2015()
    self._energy_metric = energy_metric
    self._ref_energy = None
  


  def compute_value(self, **kwargs) -> float:
    if self._ref_energy is None:
      workspace = Workspace.instance()
      import pyrosetta
      pose = pyrosetta.pose_from_pdb(workspace.target_pdb_path)
      self._ref_energy = self._energy_metric.score_fn(pose)
    other_metrics = kwargs['otherMetrics']
    e = other_metrics[self._energy_metric.column_name()] = \
      self._energy_metric(**kwargs)
    f1 = other_metrics[self._geometry_metric.column_name()] = \
      self._geometry_metric(**kwargs)
    return f1 / 4 + 1 / np.exp(abs(self._ref_energy - e) / self._weight)
  