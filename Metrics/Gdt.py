from .Metric import Metric
from typing import List
from Bio.PDB.Atom import Atom
import statistics # no se confunda con evodesign.Statistics





class Gdt(Metric):
  
  def __init__(self, cutoffs: List[float] = [ 1.0, 2.0, 4.0, 8.0 ]) -> None:
    super().__init__()
    self._cutoffs = cutoffs
  


  def __call__(self, 
               modelBackbone: List[Atom], 
               referenceBackbone: List[Atom]
               ) -> float:
    # suponemos que `modelBackbone` ya fue superpuesto contra 
    # `referenceBackbone` en un paso anterior
    distances = [ a - b for a, b in zip(modelBackbone, referenceBackbone) ]
    ratios = [ 
      sum([ d <= c for d in distances ]) / len(distances)
      for c in self._cutoffs 
    ]
    return statistics.fmean(ratios)
  