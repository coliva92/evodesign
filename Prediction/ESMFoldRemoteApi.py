from .Predictor import Predictor
import requests





class ESMFoldRemoteApi(Predictor):

  _MAX_REQUESTS = 256
  _requests_count = 0



  @classmethod
  def get_name(cls) -> str:
    return 'Predictor_ESMFold_RemoteApi'
  

  
  def predict_structure(self, 
                        sequence: str, 
                        pdbFilename: str
                        ) -> None:
    if ESMFoldRemoteApi._requests_count >= ESMFoldRemoteApi._MAX_REQUESTS:
      raise RecursionError('ESMFold API max requests reached')
    ESMFoldRemoteApi._requests_count += 1
    response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', 
                             data=sequence)
    if response.status_code == 405:
      raise RuntimeError('Forbidden')
    elif response.status_code == 408:
      raise RuntimeError('Request Timeout')
    elif response.status_code == 500:
      raise RuntimeError('Internal Server Error')
    elif response.status_code != 200:
      raise RuntimeError('Unknown HTTP error')
    with open(pdbFilename, 'wt', encoding='utf-8') as pdb_file:
        pdb_file.write(response.content.decode())
