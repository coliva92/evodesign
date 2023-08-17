from .Metric import Metric
from typing import List
from Bio.PDB.Atom import Atom
import statistics





class Gdt(Metric):
  
  def __init__(self, 
               cutoffs: List[float] = [ 1.0, 2.0, 4.0, 8.0 ],
               alphaCarbonOnly: bool = True
               ) -> None:
    super().__init__()
    self._cutoffs = cutoffs
    self._alpha_carbon_only = alphaCarbonOnly
  


  def _filter_alpha_carbon_atoms(self, backbone: List[Atom]) -> List[Atom]:
    return filter(lambda atom: atom.get_name() == 'CA', backbone)
  


  def __call__(self, 
               modelBackbone: List[Atom], 
               referenceBackbone: List[Atom]
               ) -> float:
    # suponemos que `modelBackbone` ya fue superpuesto contra 
    # `referenceBackbone` en un paso anterior
    if self._alpha_carbon_only:
      modelBackbone = self._filter_alpha_carbon_atoms(modelBackbone)
      referenceBackbone = self._filter_alpha_carbon_atoms(referenceBackbone)
    distances = [ a - b for a, b in zip(modelBackbone, referenceBackbone) ]
    return statistics.fmean([
      sum([ d <= c for d in distances ]) / len(distances)
      for c in self._cutoffs
    ])
  