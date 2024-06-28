from .Predictor import Predictor
from ..Exceptions import *
import requests
import time





class ESMFoldRemoteApi(Predictor):
  
  def predict_structure(self, 
                        sequence: str, 
                        pdbPath: str
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
    pdbPath : str
        The path to the PDB file where the predicted structure will
        be stored. 

    Raises
    ------
    Exceptions.HttpBadRequest
        Raises this exception when the API responds with HTTP error code 400.
    Exceptions.HttpForbidden
        Raises this exception when the API responds with HTTP error code 403.
    Exceptions.HttpGatewayTimeout
        Raises this exception when the API responds with HTTP error code 504.
    Exceptions.HttpInternalServerError
        Raises this exception when the API responds with HTTP error code 500.
    Exceptions.HttpUnknownError
        Raises this exception when the API responds with any other HTTP error 
        code not listed here
    """
    time.sleep(1.5)
    response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', 
                             data=sequence, 
                             timeout=30, 
                             verify=False)
    if response.status_code == 400:
      raise HttpBadRequest
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
    with open(pdbPath, 'wt', encoding='utf-8') as pdb_file:
      pdb_file.write(response.content.decode())
