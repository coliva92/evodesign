from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
import numpy as np
from typing import Optional, List


def find_atoms_in_residue(
    residue: Residue,
    atom_names: List[str],
) -> List[Optional[Atom]]:
    residue_atoms = {atom.get_name(): atom for atom in residue.get_atoms()}
    atoms = []
    for name in atom_names:
        if name in residue_atoms:
            atoms.append(residue_atoms[name])
            continue
        atoms.append(None)
    return atoms


def compute_dihedral_angle(a: Atom, b: Atom, c: Atom, d: Atom) -> float:
    # Parsons J et al. (2005) J Comput Chem 26(10):1063-1068
    u = np.cross((c.coord - b.coord), (a.coord - b.coord))
    n1 = u / np.linalg.norm(u)
    v = np.cross((d.coord - c.coord), (c.coord - b.coord))
    n2 = v / np.linalg.norm(v)
    dihedral = np.arccos(-np.dot(n1, n2))
    proyection = np.dot((c.coord - b.coord), np.cross(n1, n2))
    if proyection < 0.0:
        dihedral = -dihedral
    return dihedral
