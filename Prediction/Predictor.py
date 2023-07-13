from abc import ABC, abstractmethod
from typing import List
from Bio.PDB.Atom import Atom
from Bio.PDB import PDBParser
import evodesign.Chain as Chain
import os





class Predictor(ABC):
  """
  La representación de un algoritmo de predicción de la estructura de una 
  proteína.
  """
  
  def __init__(self) -> None:
    super().__init__()
  


  @classmethod
  @abstractmethod
  def get_name(cls) -> str:
    pass


  
  @abstractmethod
  def predict_structure(self, 
                        sequence: str, 
                        pdbFilename: str) -> None:
    """
    Predice la estructura de la secuencia de aminiácidos especificada por 
    `sequence` y escribe en resultado en un archivo PDB cuyo nombre está 
    especificado por `pdbFilename`.
    """
    pass



  def get_predicted_backbone(self, 
                             id: str,
                             sequence: str, 
                             pdbFilename: str) -> List[Atom]:
    """
    Retorna el esqueleto de la estructura predicha para la secuencia 
    especificada. 
    - `id`: el nombre que se utiliza para identificar la estructura en 
      BioPython y en el algoritmo evolutivo.
    - `sequence`: la secuencia de aminoácidos cuya estructura se desea obtener.
    - `pdbFilename`: el nombre del archivo PDB donde se guarda la estrcutura. 
      Si el archivo existe, la estructura se carga desde este archivo. En caso 
      contrario, se ejecuta el algoritmo predictor y la estructura predicha se 
      guarda en un nuevo archivo, usando este mismo nombre.
    """
    # PDBParser no puede generar una estructura directamente de los datos 
    # crudos; solo puede hacerlo desde un archivo PDB
    if not os.path.isfile(pdbFilename):
      os.makedirs(os.path.dirname(os.path.abspath(pdbFilename)), exist_ok=True)
      self.predict_structure(sequence, pdbFilename)
    parser = PDBParser()
    structure = parser.get_structure(id, pdbFilename)
    return Chain.filter_backbone_atoms_from_chain(structure)
