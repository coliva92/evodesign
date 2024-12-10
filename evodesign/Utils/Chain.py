from typing import List
from Bio.PDB.Structure import Structure
from Bio.PDB.Atom import Atom
from Bio.PDB import PDBParser
from typing import Optional
import os


BACKBONE_ATOMS = dict.fromkeys(["N", "CA", "C", "O"])


def load_structure(pdb_path: str, parser: Optional[PDBParser] = None) -> Structure:
    if parser is None:
        parser = PDBParser()
    structure_id = os.path.splitext(os.path.basename(pdb_path))[0]
    return parser.get_structure(structure_id, pdb_path)


def length(structure: Structure, model_id: int = 0, chain_id: str = "A") -> int:
    return len(list(structure[model_id][chain_id].get_residues()))


def backbone_atoms(
    structure: Structure, model_id: int = 0, chain_id: str = "A"
) -> List[Atom]:
    return [
        atom
        for atom in structure[model_id][chain_id].get_atoms()
        if atom.get_name() in BACKBONE_ATOMS
    ]


def ca_atoms(
    structure: Structure, model_id: int = 0, chain_id: str = "A"
) -> List[Atom]:
    return [
        atom
        for atom in structure[model_id][chain_id].get_atoms()
        if atom.get_name() == "CA"
    ]
