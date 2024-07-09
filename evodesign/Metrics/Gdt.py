from .Metric import Metric
from typing import List
from .Rmsd import Rmsd
import numpy as np





class Gdt(Metric):
  
  def column_name(cls) -> str:
    return 'gdt'
  

  
  def _params(self) -> dict:
    return {
      'cutoffs': self._cutoffs,
      'rmsdCalculator': self._rmsd_calc.settings()
    }
  


  def __init__(self, 
               rmsdCalculator: Rmsd,
               cutoffs: List[float] = [ 1.0, 2.0, 4.0, 8.0 ]
               ) -> None:
    """
    Computes the GDT value between the model backbone and a reference backbone
    after superimposing the former into the latter. Both backbones must contain 
    equal number of atoms.

    Parameters
    ----------
    cutoffs : List[float], optional
        The cut-off distances used for computing the GDT. Default value is 
        `[ 1.0, 2.0, 4.0, 8.0 ]`
    """
    super().__init__()
    self._cutoffs = cutoffs
    self._rmsd_calc = rmsdCalculator
  


  def compute_value(self, **kwargs) -> float:
    """
    Uses the given model and reference backbones to compute the GDT after
    superimposing the former into the latter.

    Parameters
    ----------
    model : List[Bio.PDB.Atom.Atom]
        The model backbone.
    reference : List[Bio.PDB.Atom.Atom]
        The reference backbone.
    otherMetrics : Dict[str, Union[int, float, str]], optional
        The previously computed values of other metrics. Default is an empty 
        dictionary.

    Returns
    -------
    float
        The computed GDT; it's value is between 0 and 1.
    """
    c = self._rmsd_calc.column_name()
    # check if the structures were already superimposed
    if c not in kwargs['otherMetrics']:
      kwargs['otherMetrics'][c] = self._rmsd_calc.compute_value(**kwargs)
    model, reference = kwargs['model'], kwargs['reference']
    distances = np.array([ a - b for a, b in zip(model, reference) ])
    gdt = np.mean([
      np.mean([ d <= c for d in distances ])
      for c in self._cutoffs
    ])
    return gdt
  