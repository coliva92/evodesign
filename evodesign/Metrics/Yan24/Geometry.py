from ..Metric import Metric
from ..Normalization.Reciprocal import Reciprocal
from ..Rmsd import Rmsd
from ..Gdt import Gdt
from ..ContactMapRms import ContactMapRms
from ..TMScore import TMScore
from typing import List, Optional





class Geometry(Metric):
  
  def _params(self) -> dict:
    params = super()._params()
    params['rmsd_metric'] = self._rmsd_metric.settings()
    params['gdt_metric'] = self._gdt_metric.settings()
    params['tmscore_metric'] = self._tmscore_metric.settings()
    params['cm_rms_metric'] = self._cm_rms_metric.settings()
    return params
  


  def __init__(self, 
               rmsd_metric: Optional[Rmsd] = None,
               gdt_metric: Optional[Gdt] = None,
               tmscore_metric: Optional[TMScore] = None,
               cm_rms_metric: Optional[ContactMapRms] = None, 
               column: Optional[str] = None
               ) -> None:
    super().__init__(column)
    self._rmsd_metric = rmsd_metric
    if self._rmsd_metric is None:
      self._rmsd_metric = Rmsd()
    self._gdt_metric = gdt_metric
    if self._gdt_metric is None:
      self._gdt_metric = Gdt(rmsd_metric=self._rmsd_metric)
    self._tmscore_metric = tmscore_metric
    if self._tmscore_metric is None:
      self._tmscore_metric = TMScore(rmsd_metric=self._rmsd_metric)
    self._cm_rms_metric = cm_rms_metric
    if self._cm_rms_metric is None:
      self._cm_rms_metric = ContactMapRms()
  


  def compute_value(self, **kwargs) -> float:
    other_metrics = kwargs['otherMetrics']
    a = self._rmsd_metric(**kwargs)
    b = other_metrics[self._gdt_metric.column_name()] = self._gdt_metric(**kwargs)
    c = self._cm_rms_metric(**kwargs)
    d = other_metrics[self._tmscore_metric.column_name()] = \
      self._tmscore_metric(**kwargs)
    return a + b + c + d
  