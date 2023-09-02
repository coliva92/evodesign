from .Metric import Metric
from typing import List, Optional
from Bio.PDB.Atom import Atom
import evodesign.Chain as Chain
import statistics





class Gdt(Metric):
  
  def __init__(self, 
               cutoffs: List[float] = [ 1.0, 2.0, 4.0, 8.0 ],
               carbonAlphaOnly: bool = False
               ) -> None:
    super().__init__()
    self._cutoffs = cutoffs
    self._carbon_alpha_only = carbonAlphaOnly
  


  def __call__(self, 
               modelBackbone: List[Atom], 
               referenceBackbone: List[Atom],
               _: Optional[str] = None
               ) -> float:
    # suponemos que `modelBackbone` ya fue superpuesto contra 
    # `referenceBackbone` en un paso anterior
    if self._carbon_alpha_only:
      modelBackbone = Chain.filter_alpha_carbons_in_backbone(modelBackbone)
      referenceBackbone = Chain.filter_alpha_carbons_in_backbone(
        referenceBackbone)
    distances = [ a - b for a, b in zip(modelBackbone, referenceBackbone) ]
    # tradicionalmente, el GDT se reporta en centésimas
    return 100 * statistics.fmean([
      sum([ d <= c for d in distances ]) / len(distances)
      for c in self._cutoffs
    ])
  