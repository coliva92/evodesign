from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure


@dataclass
class Chain:
    sequence: str | None = None
    sequence_numpy: npt.NDArray[np.float64] | None = None
    structure: Structure | None = None
    model_id: int | None = None
    chain_id: str | None = None
    pdb_path: str | None = None
    residues: list[Residue] | None = None
    backbone_atoms: list[Atom] | None = None
    ca_atoms: list[Atom] | None = None
