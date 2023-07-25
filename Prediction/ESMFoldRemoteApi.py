from .Predictor import Predictor
import requests





class ESMFoldRemoteApi(Predictor):
  """
  El algoritmo ESMFold para predecir la estructura de una proteína.
  """

  _MAX_REQUESTS = 256
  _requests_count = 0



  def __init__(self) -> None:
    super().__init__()
  


  @classmethod
  def get_name(cls) -> str:
    return 'Predictor_ESMFold_RemoteApi'
  

  
  def predict_structure(self, 
                        sequence: str, 
                        pdbFilename: str) -> None:
    """
    Predice la estructura de la secuencia de aminiácidos especificada por 
    `sequence` y escribe en resultado en un archivo PDB cuyo nombre está 
    especificado por `pdbFilename`.
    """
    if ESMFoldRemoteApi._requests_count >= ESMFoldRemoteApi._MAX_REQUESTS:
      raise RecursionError('ESMFold API max requests reached')
    ESMFoldRemoteApi._requests_count += 1
    response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', 
                             data=sequence)
    if response.status_code == 405:
      raise RuntimeError('Forbidden')
    if response.status_code == 408:
      raise RuntimeError('Request Timeout')
    if response.status_code == 500:
      raise RuntimeError('Internal Server Error')
    if response.status_code != 200:
      raise RuntimeError('Unknown HTTP error')
    with open(pdbFilename, 'wt', encoding='utf-8') as the_file:
        the_file.write(response.content.decode())
