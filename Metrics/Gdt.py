from .Metric import Metric
from typing import List, Optional
from Bio.PDB.Atom import Atom
import statistics





class Gdt(Metric):
  
  def __init__(self, 
               cutoffs: List[float] = [ 1.0, 2.0, 4.0, 8.0 ]
               ) -> None:
    super().__init__()
    self._cutoffs = cutoffs
  


  def __call__(self, 
               modelBackbone: List[Atom], 
               referenceBackbone: List[Atom],
               _: Optional[str] = None
               ) -> float:
    # suponemos que `modelBackbone` ya fue superpuesto contra 
    # `referenceBackbone` en un paso anterior
    distances = [ a - b for a, b in zip(modelBackbone, referenceBackbone) ]
    # tradicionalmente, el GDT se reporta en centésimas
    return statistics.fmean([
      sum([ d <= c for d in distances ]) / len(distances)
      for c in self._cutoffs
    ])
  