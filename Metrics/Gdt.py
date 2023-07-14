from .Metric import Metric
from typing import List, Optional
from Bio.PDB.Atom import Atom
from Bio.PDB import Superimposer





class Gdt(Metric):
  """
  El RMSD obtenido de una superposición de dos estructuras.
  """
  
  def __init__(self, cutoffs: List[float] = [ 1.0, 2.0, 4.0, 8.0 ]) -> None:
    super().__init__()
    self.cutoffs = cutoffs
  


  def compute(self, 
              model_backbone: List[Atom],
              reference_backbone: Optional[List[Atom]]) -> float:
    # suponemos que `model_backbone` ya fue superpuesto contra 
    # `reference_backbone` en un paso anterior
    n = len(reference_backbone)
    distances = [
      model_backbone[i] - reference_backbone[i] for i in range(n)
    ]
    ratios_sum = 0.0
    for cutoff in self.cutoffs:
      tally = 0
      for d in distances:
        tally += d <= cutoff
      ratios_sum += tally / n
    return ratios_sum / 4.0
  