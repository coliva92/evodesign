from Bio.PDB.Atom import Atom
from .Metric import Metric
from typing import List, Optional, Tuple
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
    
    def preserved(i: int, j: int, cutoff: float) -> bool:
      a = modelBackbone[i]
      b = modelBackbone[j]
      x = referenceBackbone[i]
      y = referenceBackbone[j]
      return abs((x - y) - (a - b)) <= cutoff
    
    reference_set = {
      frozenset([ i, j ])
      for i, a in enumerate(referenceBackbone)
      for j, b in enumerate(referenceBackbone)
      if a - b > 0 and a - b <= self._radius
    }
    model_distances = [
        ( i, j )
        for i, a in enumerate(modelBackbone)
        for j, b in enumerate(modelBackbone)
        if a - b > 0 and a - b <= self._radius and { i, j } in reference_set
    ]
    
    n = len(reference_set)
    return statistics.fmean([
      sum([ preserved(x[0], x[1], c) for x in model_distances ]) / n
      for c in self._cutoffs
    ])
