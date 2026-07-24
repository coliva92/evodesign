from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
import numpy as np
import numpy.typing as npt
import itertools
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


def compute_distance_map(atoms: List[Atom]) -> npt.NDArray[np.float64]:
    return np.array([a - b for (a, b) in itertools.combinations(atoms, 2)])


def compute_contacts_map(atoms: List[Atom], 
                         distance_threshold: float = 8.0,
                         min_separation: int = 6,
                         ) -> npt.NDArray[np.int64]:
    L = len(atoms)
    if L == 0:
        raise ValueError
    distances = compute_distance_map(atoms)
    row_indices, column_indices = np.triu_indices(L, k=1)
    separation_mask = (column_indices - row_indices) >= min_separation
    contacts_mask = distances <= distance_threshold
    contacts = separation_mask & contacts_mask
    contacts_map = np.zeros((L, L), dtype=np.int64)
    rows = row_indices[contacts]
    columns = column_indices[contacts]
    contacts_map[rows, columns] = 1
    contacts_map[columns, rows] = 1 
    return contacts_map
