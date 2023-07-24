from .Predictor import Predictor
import requests
import time





class ESMFold(Predictor):
  """
  El algoritmo ESMFold para predecir la estructura de una proteína.
  """
  
  forbiddens_count = 0
  server_errors_count = 0



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
    if response.status_code == 405:
      ESMFold.forbiddens_count += 1
      raise RuntimeError('Forbidden')
    if response.status_code == 500:
      ESMFold.server_errors_count += 1
      raise RuntimeError('Internal server error')
    if response.status_code != 200:
      raise RuntimeError('Unknown HTTP error')
    ESMFold.forbiddens_count, ESMFold.server_errors_count = 0, 0
    with open(pdbFilename, 'wt', encoding='utf-8') as file:
        file.write(response.content.decode())



def handle_api_errors(e: RuntimeError) -> None:
  if e == 'Forbidden': 
    sleep_time = 620 if ESMFold.forbiddens_count == 1 else 3620
    time.sleep(sleep_time)
  if e == 'Internal server error' and \
      ESMFold.server_errors_count >= 5:
    time.sleep(300)
