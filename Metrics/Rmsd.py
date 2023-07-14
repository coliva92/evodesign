from .Metric import Metric
from typing import List
from Bio.PDB.Atom import Atom
from Bio.PDB import Superimposer





class Rmsd(Metric):
  """
  El RMSD obtenido de una superposición de dos estructuras.
  """
  
  def __init__(self) -> None:
    super().__init__()
  


  def compute(self, 
              model_backbone: List[Atom],
              reference_backbone: List[Atom]) -> float:
    """
    Superpone el esqueleto especificado por `backbone` sobre el esqueleto de 
    referencia y retorna el RMSD resultante.
    """
    superimposer = Superimposer()
    superimposer.set_atoms(reference_backbone, model_backbone)
    superimposer.apply(model_backbone)
    return superimposer.rms
  