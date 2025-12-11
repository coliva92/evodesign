from .Predictor import Predictor
from .DirectoryManager import DirectoryManager
from Bio.PDB import PDBIO, Structure, Model, Chain, Residue, Atom
from Bio.Data import IUPACData
from typing import List
import os


class Null(Predictor):

    _io = PDBIO()
    _amino_acids = {
        k: v.capitalize() for k, v in IUPACData.protein_letters_1to3.items()
    }

    def predict_single_pdb_file(
        self,
        sequence: str,
        protein_name: str,
        directory: DirectoryManager,
    ) -> None:
        # Create a dummy PDB to pass the amino acid sequence to the next component
        structure = Structure.Structure(protein_name)
        model = Model.Model(0)
        chain = Chain.Chain("A")
        for j, aa in enumerate(sequence, start=1):
            resname = self._amino_acids.get(aa, "UNK")
            res_id = (" ", j, " ")
            residue = Residue.Residue(res_id, resname.upper(), "")
            residue.add(
                Atom.Atom("N", [j * 1.5, 0.0, 0.0], 1.0, 1.0, " ", "N", j, element="N")
            )
            residue.add(
                Atom.Atom(
                    "CA", [j * 1.5, 1.5, 0.0], 1.0, 1.0, " ", "CA", j, element="C"
                )
            )
            residue.add(
                Atom.Atom("C", [j * 1.5, 3.0, 0.0], 1.0, 1.0, " ", "C", j, element="C")
            )
            residue.add(
                Atom.Atom("O", [j * 1.5, 4.0, 0.0], 1.0, 1.0, " ", "O", j, element="O")
            )
            chain.add(residue)
        model.add(chain)
        structure.add(model)
        self._io.set_structure(structure)
        output_path = os.path.join(directory.prediction_pdbs_dir, f"{protein_name}.pdb")
        self._io.save(output_path)

    def do(
        self,
        sequences: List[str],
        directory: DirectoryManager,
    ) -> None:
        os.makedirs(directory.prediction_pdbs_dir, exist_ok=True)
        directory.empty_folders_content()
        for i, sequence in enumerate(sequences):
            protein_name = f"{directory.prefix}_{i}"
            self.predict_single_pdb_file(sequence, protein_name, directory)
