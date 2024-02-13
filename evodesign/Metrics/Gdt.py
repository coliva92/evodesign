from .Metric import Metric
from typing import List
import numpy as np





class Gdt(Metric):
  
  @classmethod
  def column_name(cls) -> str:
    return 'gdt'
  

  
  def __init__(self, 
               cutoffs: List[float] = [ 1.0, 2.0, 4.0, 8.0 ]
               ) -> None:
    """
    This class computes the GDT value between two backbones of equal length.

    Parameters
    ----------
    cutoffs : List[float], optional
        The cut-off distances used for computing the GDT. Default value is 
        [ 1.0, 2.0, 4.0, 8.0 ]
    """
    super().__init__()
    self.cutoffs = cutoffs
  


  def __call__(self, **kwargs) -> float:
    """
    Computes the RMSD between the atom coordinates of a model backbone and
    a reference backbone. Both backbones must have equal number of atoms.

    Parameters
    ----------
    model : List[Bio.PDB.Atom.Atom]
        The model backbone.
    reference : List[Bio.PDB.Atom.Atom]
        The reference backbone.

    Returns
    -------
    float
        The computed GDT; it's value is between 0 and 1.
    """
    # we assume that the backbones are already superimposed
    model, reference = kwargs['model'], kwargs['reference']
    distances = np.array([ a - b for a, b in zip(model, reference) ])
    return np.mean([
      np.mean([ d <= c for d in distances ])
      for c in self.cutoffs
    ])
  