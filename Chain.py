"""Colección de funciones auxiliares para trabajar con una instancia de 
`Bio.PDB.Chain` de BioPython, la cual representa una cadena de aminoácidos.
"""
from typing import List
from Bio.PDB.Structure import Structure
from Bio.PDB.Atom import Atom
from Bio.PDB import PDBParser
import os





"""Las letras válidas que representan los átomos del esqueleto polipeptídico.
"""
BACKBONE_ATOMS = [ 'N', 'CA', 'C' ]



def load_structure_from_pdb_file(filename: str) -> Structure:
  """Carga una proteína desde el archivo PDB especificado por `filename`.
  """
  id = os.path.splitext(os.path.basename(filename))[0]
  parser = PDBParser()
  return parser.get_structure(id, filename)



def count_residues_from_chain(structure: Structure, 
                              modelId: int = 0, 
                              chainId: str = 'A'
                              ) -> int:
  """Retorna el número de residuos de una cadena en la proteína especificada.
  Se requiere conocer el identificador del modelo y la cadena para ejecutar 
  esta función (ver la sección "The Structure Object" de [esta documentación](https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ)).
  - `structure`: la proteína cuyo número de residuos se desea contar.
  - `modelId`: el identificador del modelo que contiene la cadena cuyos 
    residuos se desean contar.
  - `chainId`: el identificador de la cadena cuyos residuos se desean contar.
  """
  return len(list(structure[modelId][chainId].get_residues()))



def filter_backbone_atoms_from_chain(structure: Structure,
                                     modelId: int = 0, 
                                     chainId: str = 'A'
                                     ) -> List[Atom]:
  """Retorna un arreglo que contiene únicamente los átomos del esqueleto 
  polipeptídico de la proteína especificada. Se requiere conocer el 
  identificador del modelo y la cadena para ejecutar 
  esta función (ver la sección "The Structure Object" de [esta documentación](https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ)).
  - `structure`: la proteína cuyos átomos van a ser filtrados.
  - `modelId`: el identificador del modelo que contiene la cadena cuyos 
    átomos van a ser filtrados.
  - `chainId`: el identificador de la cadena cuyos átomos van a ser filtrados.
  """
  backbone = []
  for atom in structure[modelId][chainId].get_atoms():
    if atom.get_name() in BACKBONE_ATOMS:
      backbone.append(atom)
  return backbone
