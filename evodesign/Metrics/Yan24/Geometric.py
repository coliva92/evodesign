from ..Metric import Metric
from ..Normalization.Reciprocal import Reciprocal
from ..Rmsd import Rmsd
from ..Gdt import Gdt
from ..ContactMapRms import ContactMapRms
from ..TMScore import TMScore
from typing import List





class Geometric(Metric):

  @classmethod
  def _class_name(cls) -> str:
    return 'Metrics.Yan24.Geometric'
  


  def column_name(self) -> str:
    return 'yan24_geometric'
  


  def _params(self) -> dict:
    return { 'cutoffs': self._cutoffs }
  


  def __init__(self, 
               cutoffs: List[float] = [ 1.0, 2.0, 4.0, 8.0 ]
               ) -> None:
    super().__init__()
    temp = Rmsd()
    self._rmsd_calc = Reciprocal(temp)
    self._gdt_calc = Gdt(temp, cutoffs)
    self._contact_map_calc = Reciprocal(ContactMapRms())
    self._tm_score_calc = TMScore(temp)
    self._cutoffs = cutoffs
  


  def compute_value(self, **kwargs) -> float:
    other_metrics = kwargs['otherMetrics']
    a = self._rmsd_calc(**kwargs)
    b = other_metrics[self._gdt_calc.column_name()] = self._gdt_calc(**kwargs)
    c = self._contact_map_calc(**kwargs)
    d = other_metrics[self._tm_score_calc.column_name()] = \
      self._tm_score_calc(**kwargs)
    return a + b + c + d
  