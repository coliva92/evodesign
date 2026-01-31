from .AminoAcids import to_numpy, to_str
from typing import List, Optional
from .Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.PDB.Atom import Atom
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import three_to_one, is_aa
from Bio.PDB.Residue import Residue
import numpy as np
import numpy.typing as npt
import os


BACKBONE_ATOMS = dict.fromkeys(["N", "CA", "C", "O"])


class ChainFactory:

    _parser = None

    @classmethod
    def create_from_pdb(
        cls,
        pdb_path: str,
        model_id: int = 0,
        chain_id: Optional[str] = None,
    ) -> Chain:
        if cls._parser is None:
            cls._parser = PDBParser(QUIET=True)
        pdb_path = os.path.abspath(pdb_path)
        structure_id = os.path.splitext(os.path.basename(pdb_path))[0]
        structure = cls._parser.get_structure(structure_id, pdb_path)
        if chain_id is None:
            for chain_id in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
                if chain_id in structure[model_id]:
                    break
        sequence = cls.get_sequence(structure, model_id, chain_id)
        sequence_numpy = cls.sequence_str_to_numpy(sequence)
        backbone_atoms = cls.get_backbone_atoms(structure, model_id, chain_id)
        ca_atoms = cls.get_ca_atoms(structure, model_id, chain_id)
        residues = cls.get_residues(structure, model_id, chain_id)
        return Chain(
            sequence,
            sequence_numpy,
            structure,
            model_id,
            chain_id,
            pdb_path,
            residues,
            backbone_atoms,
            ca_atoms,
        )

    @classmethod
    def create_from_numpy(cls, sequence_numpy: npt.NDArray[np.int64]) -> Chain:
        sequence = cls.sequence_numpy_to_str(sequence_numpy)
        return Chain(sequence=sequence, sequence_numpy=sequence_numpy)

    @classmethod
    def get_sequence(cls, structure: Structure, model_id: int, chain_id: str) -> str:
        chain = structure[model_id][chain_id]
        return "".join(three_to_one(r.get_resname()) for r in chain if is_aa(r))

    @classmethod
    def get_backbone_atoms(
        cls, structure: Structure, model_id: int, chain_id: str
    ) -> List[Atom]:
        return [
            atom
            for atom in structure[model_id][chain_id].get_atoms()
            if atom.get_name() in BACKBONE_ATOMS
        ]

    @classmethod
    def get_ca_atoms(
        cls, structure: Structure, model_id: int, chain_id: str
    ) -> List[Atom]:
        return [
            atom
            for atom in structure[model_id][chain_id].get_atoms()
            if atom.get_name() == "CA"
        ]

    @classmethod
    def get_residues(
        cls, structure: Structure, model_id: int, chain_id: str
    ) -> List[Residue]:
        return [res for res in structure[model_id][chain_id].get_residues()]

    @classmethod
    def sequence_str_to_numpy(cls, sequence: str) -> npt.NDArray[np.int64]:
        return to_numpy(sequence)

    @classmethod
    def sequence_numpy_to_str(cls, sequence_numpy: npt.NDArray[np.int64]) -> str:
        return to_str(sequence_numpy)
