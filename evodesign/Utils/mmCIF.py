from Bio.PDB import MMCIFParser, PDBIO


def convert_cif_to_pdb(
    cif_path: str,
    structure_id: str = "A",
    parser: MMCIFParser = MMCIFParser(),
    io: PDBIO = PDBIO(),
) -> str:
    structure = parser.get_structure(structure_id, cif_path)
    pdb_path = cif_path.replace(".cif", ".pdb")
    io.set_structure(structure)
    io.save(pdb_path)
    return pdb_path
