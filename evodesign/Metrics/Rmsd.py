from .Metric import Metric
from Bio.PDB import Superimposer





class Rmsd(Metric):

  _superimposer = None
  


  def column_name(cls) -> str:
    return 'rmsd'

  

  def compute_value(self, **kwargs) -> float:
    """
    Computes the root-mean-squared of the distance between the atom coordinates 
    of a model backbone and a reference backbone after superimposing the former 
    into the latter. Both backbones must contain equal number of atoms.

    Parameters
    ----------
    model : List[Bio.PDB.Atom.Atom]
        The model backbone. In this case, the model backbone will be moved 
        towards the reference backbone during the superposition.
    reference : List[Bio.PDB.Atom.Atom]
        The reference backbone. In this case, the reference backbone remains
        fixed in position during the superposition.

    Returns
    -------
    float
        The computed RMSD measured in Angstroms.
    """
    if not self._superimposer:
      self._superimposer = Superimposer()
    self._superimposer.set_atoms(kwargs['reference'], kwargs['model'])
    self._superimposer.apply(kwargs['model'])
    return self._superimposer.rms
  