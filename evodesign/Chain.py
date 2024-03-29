from typing import List
from Bio.PDB.Structure import Structure
from Bio.PDB.Atom import Atom
from Bio.PDB import PDBParser
import os





class Chain:

  BACKBONE_ATOMS = dict.fromkeys([ 'N', 'CA', 'C', 'O' ])
  _pdb_parser = None



  @classmethod
  def length(cls,
             structure: Structure, 
             modelId: int = 0, 
             chainId: str = 'A'
             ) -> int:
    """
    Get the number of residues in a given polypeptide chain.

    Parameters
    ----------
    structure : Bio.PDB.Structure.Structure
        The BioPython Structure containing the desired polypeptide. Structure 
        objects in BioPython may contain multiple polypeptide chains, organized
        in a hierarchical representation. However, this function can only work
        with a single chain at a time, which must be specified by means of two
        IDs--one for the model and one for the actual chain (see "The Structure
        Object" section provided here: https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ).
    modelId : int, optional
        The model ID for finding the desired chain. The default is 0.
    chainId : str, optional
        The chain ID of the desired chain. The default is 'A'.

    Returns
    -------
    int
        The number of residues in the specified chain.
    """
    return len(list(structure[modelId][chainId].get_residues()))



  @classmethod
  def backbone_atoms(cls,
                     structure: Structure,
                     modelId: int = 0, 
                     chainId: str = 'A'
                     ) -> List[Atom]:
    """
    Get the backbone atoms of a given polypeptide chain.

    Parameters
    ----------
    structure : Bio.PDB.Structure.Structure
        The BioPython Structure containing the desired polypeptide. Structure 
        objects in BioPython may contain multiple polypeptide chains, organized
        in a hierarchical representation. However, this function can only work
        with a single chain at a time, which must be specified by means of two
        IDs--one for the model and one for the actual chain (see "The Structure
        Object" section provided here: https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ).
    modelId : int, optional
        The model ID for finding the desired chain. The default is 0.
    chainId : str, optional
        The chain ID of the desired chain. The default is 'A'.

    Returns
    -------
    List[Atom]
        The backbone atoms of the specified chain.
    """
    return [
      atom
      for atom in structure[modelId][chainId].get_atoms()
      if atom.get_name() in cls.BACKBONE_ATOMS
    ]
  


  @classmethod
  def load_structure(cls, pdbPath: str) -> Structure:
    """
    Loads a BioPython Structure instance from the specified PDB file.

    Parameters
    ----------
    pdbPath : str
        The path to the PDB file to be opened.

    Returns
    -------
    Bio.PDB.Structure.Structure
        The loaded Structure instance.
    """
    structure_id = os.path.splitext(os.path.basename(pdbPath))[0]
    if not cls._pdb_parser:
      cls._pdb_parser = PDBParser()
    return cls._pdb_parser.get_structure(structure_id, pdbPath)
