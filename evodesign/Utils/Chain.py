from typing import List, Optional
from Bio.PDB.Structure import Structure
from Bio.PDB.Atom import Atom
from Bio.PDB import PDBParser, PPBuilder
import numpy as np
import numpy.typing as npt
import os


def numpy_sequence_to_str(sequence_numpy: npt.NDArray[np.int64]):
    return "".join(Chain.AMINO_ACIDS[i] for i in sequence_numpy)


class Chain:

    _ppb = None
    BACKBONE_ATOMS = dict.fromkeys(["N", "CA", "C", "O"])
    AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
    MAP_AMINO_ACID_TO_INT = {aa: i for i, aa in enumerate(AMINO_ACIDS)}

    def __init__(
        self,
        structure: Structure,
        model_id: int = 0,
        chain_id: str = "A",
        pdb_path: Optional[str] = None
    ):
        self.structure = structure
        self.model_id = model_id
        self.chain_id = chain_id
        self.pdb_path = pdb_path
        self.sequence = self._get_sequence()
        self.sequence_numpy = np.array(
            [self.MAP_AMINO_ACID_TO_INT[aa] for aa in self.sequence], dtype=np.int64
        )
        self.backbone_atoms = self._get_backbone_atoms()
        self.ca_atoms = self._get_ca_atoms()

    def _get_sequence(self) -> str:
        if self._ppb is None:
            self._ppb = PPBuilder()
        chain = self.structure[self.model_id][self.chain_id]
        return "".join(str(pp.get_sequence()) for pp in self._ppb.build_peptides(chain))

    def _get_backbone_atoms(self) -> List[Atom]:
        return [
            atom
            for atom in self.structure[self.model_id][self.chain_id].get_atoms()
            if atom.get_name() in self.BACKBONE_ATOMS
        ]

    def _get_ca_atoms(self) -> List[Atom]:
        return [
            atom
            for atom in self.structure[self.model_id][self.chain_id].get_atoms()
            if atom.get_name() == "CA"
        ]


class ChainFactory:

    _parser = None

    @classmethod
    def create(
        cls,
        pdb_path: str,
        model_id: int = 0,
        chain_id: Optional[str] = None,
    ) -> Chain:
        if cls._parser is None:
            cls._parser = PDBParser(QUIET=True)
        structure_id = os.path.splitext(os.path.basename(pdb_path))[0]
        structure = cls._parser.get_structure(structure_id, pdb_path)
        if chain_id is None:
            for chain_id in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
                if chain_id in structure[model_id]:
                    break
        return Chain(structure, model_id, chain_id, pdb_path)
