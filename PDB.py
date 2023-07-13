from Bio.PDB import PDBParser
from Bio.PDB.Structure import Structure
import os





"""
Colección de funciones auxiliares para leer un archivo de PDB (u otros 
formatos) y cargar su contenido a una instancia `Structure` de biopython.
"""

def load_structure_from_pdb_file(filename: str) -> Structure:
  """
  Carga una proteína desde el archivo PDB especificado.
    - `filename`: el nombre del archivo PDB cuyo contenido será cargado.
    - `id`: el nombre que se utilizará para identificar internamente la 
    proteína cargada.
  """
  id = os.path.splitext(os.path.basename(filename))[0]
  parser = PDBParser()
  return parser.get_structure(id, filename)
