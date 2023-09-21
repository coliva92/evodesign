from Bio.PDB.Atom import Atom
from .Metric import Metric
from typing import List, Optional
import statistics





class Lddt(Metric):

  def __init__(self,
               cutoffs: List[float] = [ 0.5, 1.0, 2.0, 4.0 ],
               inclusionRadius: float = 15,
               ) -> None:
    super().__init__()
    self._radius = inclusionRadius
    self._cutoffs = cutoffs



  def __call__(self, 
               modelBackbone: List[Atom], 
               referenceBackbone: List[Atom], 
               _: Optional[str] = None) -> float:
    reference_count = len(list(filter(lambda d: d > 0 and d <= self._radius, 
                                      [
                                        a - b 
                                        for a in referenceBackbone 
                                        for b in referenceBackbone 
                                      ])))
    model_distances = list(filter(lambda d: d > 0 and d <= self._radius,
                                  [ 
                                    a - b 
                                    for a in modelBackbone 
                                    for b in modelBackbone 
                                  ]))
    return statistics.fmean([ 
      len(list(filter(lambda d: d <= c, model_distances))) / reference_count
      for c in self._cutoffs 
    ])
