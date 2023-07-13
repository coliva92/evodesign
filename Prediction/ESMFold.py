from .Predictor import Predictor
import requests
import time
import random





class ESMFold(Predictor):
  """
  El algoritmo ESMFold para predecir la estructura de una proteína.
  """
  
  def __init__(self) -> None:
    super().__init__()
  


  @classmethod
  def get_name(cls) -> str:
    return 'Predictor_ESMFold'
  

  
  def predict_structure(self, 
                        sequence: str, 
                        pdbFilename: str) -> None:
    """
    Predice la estructura de la secuencia de aminiácidos especificada por 
    `sequence` y escribe en resultado en un archivo PDB cuyo nombre está 
    especificado por `pdbFilename`.
    """
    response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', 
                             data=sequence)
    if response.status_code != 200:
      raise RuntimeError('a call to the ESMFold remote API failed.')
    with open(pdbFilename, 'wt', encoding='utf-8') as file:
        file.write(response.content.decode())
