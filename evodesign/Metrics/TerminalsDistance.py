from . import Metric
from typing import List, Optional
from Bio.PDB.Atom import Atom





class TerminalsDistance(Metric):

  def __call__(self, 
               modelBackbone: List[Atom], 
               referenceBackbone: List[Atom],
               _: Optional[str] = None
               ) -> float:
    # calcular la distancia entre N del primer residuo y C del último
    return modelBackbone[0] - modelBackbone[-2]
