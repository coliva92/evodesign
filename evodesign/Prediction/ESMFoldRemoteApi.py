from .Predictor import Predictor
from ..Exceptions import (HttpForbidden, 
                        HttpInternalServerError, 
                        HttpGatewayTimeout,
                        HttpUnknownError)
import requests
import time





class ESMFoldRemoteApi(Predictor):

  @classmethod
  def name(cls) -> str:
    return 'Predictor_ESMFold_RemoteApi'
  

  
  def predict_structure(self, 
                        sequence: str, 
                        pdbFilename: str
                        ) -> None:
    """
    Predicts the 3D structure of a given amino acid sequence using ESMFold's
    remote API.

    Parameters
    ----------
    sequence : str
        The amino acid sequence which structure will be predicted. Each residue
        must be represented with a single letter corresponding to one of the
        20 essential amino acids.
    pdbFilename : str
        The path and name of the PDB file where the predicted structure will
        be stored. 

    Raises
    ------
    HttpForbidden
        Raises this exception when the API responds with HTTP error code 403.
    HttpGatewayTimeout
        Raises this exception when the API responds with HTTP error code 504.
    HttpInternalServerError
        Raises this exception when the API responds with HTTP error code 500.
    HttpUnknownError
        Raises this exception when the API responds with any other HTTP error 
        code.
    """
    time.sleep(1.5)
    response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', 
                             data=sequence, 
                             timeout=30, 
                             verify=False)
    if response.status_code == 403:
      raise HttpForbidden
    elif response.status_code == 504:
      raise HttpGatewayTimeout
    elif response.status_code == 500:
      raise HttpInternalServerError
    elif response.status_code != 200:
      print(response.status_code)
      print(response.content.decode())
      raise HttpUnknownError
    with open(pdbFilename, 'wt', encoding='utf-8') as pdb_file:
      pdb_file.write(response.content.decode())