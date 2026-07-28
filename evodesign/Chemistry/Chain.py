from dataclasses import dataclass
from typing import List, Optional

import numpy as np
import numpy.typing as npt
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure


@dataclass
class Chain:

    sequence: Optional[str] = None
    sequence_numpy: Optional[npt.NDArray[np.float64]] = None
    structure: Optional[Structure] = None
    model_id: Optional[int] = None
    chain_id: Optional[str] = None
    pdb_path: Optional[str] = None
    residues: Optional[List[Residue]] = None
    backbone_atoms: Optional[List[Atom]] = None
    ca_atoms: Optional[List[Atom]] = None
