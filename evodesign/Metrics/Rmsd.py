from .Metric import Metric
from typing import List
from Bio.PDB.Atom import Atom
from Bio.PDB import Superimposer





class Rmsd(Metric):

  _superimposer = None

  

  def __call__(self, 
               sequence: str,
               model: List[Atom], 
               reference: List[Atom]
               ) -> float:
    """
    Computes the RMSD between the atom coordinates of a model backbone and
    a reference backbone after superimposing the former into the latter. Both
    backbones must contain equal number of atoms.

    Parameters
    ----------
    sequence : str
        Unused.
    model : List[Atom]
        The model backbone. In this case, the model backbone will be moved 
        towards the reference backbone during the superposition.
    reference : List[Atom]
        The reference backbone. In this case, the reference backbone remains
        fixed in position during the superposition.

    Returns
    -------
    float
        The computed RMSD measured in Angstroms.
    """
    if not self._superimposer:
      self._superimposer = Superimposer()
    self._superimposer.set_atoms(reference, model)
    self._superimposer.apply(model)
    return self._superimposer.rms
  