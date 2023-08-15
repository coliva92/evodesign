from .Predictor import Predictor
from ..Exceptions import (HttpForbidden, 
                          HttpInternalServerError, 
                          HttpGatewayTimeout,
                          HttpUnknownError)
import requests





class ESMFoldRemoteApi(Predictor):

  @classmethod
  def get_name(cls) -> str:
    return 'Predictor_ESMFold_RemoteApi'
  

  
  def predict_structure(self, 
                        sequence: str, 
                        pdbFilename: str
                        ) -> None:
    response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', 
                             data=sequence, 
                             timeout=30)
    # print(response.status_code)
    if response.status_code == 403:
      raise HttpForbidden
    elif response.status_code == 504 or requests.exceptions.ConnectTimeout:
      raise HttpGatewayTimeout
    elif response.status_code == 500:
      raise HttpInternalServerError
    elif response.status_code != 200:
      print(response.status_code)
      print(response.content.decode())
      raise HttpUnknownError
    with open(pdbFilename, 'wt', encoding='utf-8') as pdb_file:
      pdb_file.write(response.content.decode())
