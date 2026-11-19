import os

from Bio.PDB import PDBIO, MMCIFParser


def convert_cif_to_pdb(
    cif_path: str,
    structure_id: str = "A",
    parser: MMCIFParser | None = None,
    io: PDBIO | None = None,
) -> str:

    if parser is None:
        parser = MMCIFParser()

    if io is None:
        io = PDBIO()

    structure = parser.get_structure(structure_id, cif_path)
    pdb_path = os.path.splitext(cif_path)[0] + ".pdb"
    io.set_structure(structure)
    io.save(pdb_path)
    return pdb_path
