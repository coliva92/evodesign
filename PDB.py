"""Colección de funciones auxiliares para leer un archivo de PDB y cargar su 
contenido a una instancia `Bio.PDB.Structure` de BioPython.
"""
from Bio.PDB import PDBParser
from Bio.PDB.Structure import Structure
import os





def load_structure_from_pdb_file(filename: str) -> Structure:
  """Carga una proteína desde el archivo PDB especificado por `filename`.
  """
  id = os.path.splitext(os.path.basename(filename))[0]
  parser = PDBParser()
  return parser.get_structure(id, filename)
