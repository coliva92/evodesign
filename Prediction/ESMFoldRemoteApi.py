from .Predictor import Predictor
from ..Exceptions import (HttpForbidden, 
                      HttpInternalServerError, 
                      HttpGatewayTimeout, 
                      RemoteApiRequestsExceeded)
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
      raise RemoteApiRequestsExceeded
    ESMFoldRemoteApi._requests_count += 1
    response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', 
                             data=sequence)
    if response.status_code == 403:
      raise HttpForbidden
    elif response.status_code == 504:
      raise HttpGatewayTimeout
    elif response.status_code == 500:
      raise HttpInternalServerError
    elif response.status_code != 200:
      print(response.status_code)
      print(response.content.decode())
      raise RuntimeError
    with open(pdbFilename, 'wt', encoding='utf-8') as pdb_file:
      pdb_file.write(response.content.decode())
