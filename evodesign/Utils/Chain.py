from dataclasses import dataclass
from typing import List, Optional
from Bio.PDB.Structure import Structure
from Bio.PDB.Atom import Atom
import numpy as np
import numpy.typing as npt


@dataclass
class Chain:

    sequence: Optional[str] = None
    sequence_numpy: Optional[npt.NDArray[np.float64]] = None
    structure: Optional[Structure] = None
    model_id: Optional[int] = None
    chain_id: Optional[str] = None
    pdb_path: Optional[str] = None
    backbone_atoms: Optional[List[Atom]] = None
    ca_atoms: Optional[List[Atom]] = None
