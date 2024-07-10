from .Predictor import Predictor
from ..Exceptions import *
import requests
import time
from typing import Optional





class ESMFoldRemoteApi(Predictor):

    def _params(self) -> dict:
        return { 
            'url': self._url,
            'request_is_json': self._request_is_json,
            'response_is_json': self._response_is_json
        }
    


    def __init__(self, 
                 url: str = 'https://api.esmatlas.com/foldSequence/v1/pdb/',
                 request_is_json: bool = False,
                 response_is_json: bool = False
                 ) -> None:
        """
        Interface for interacting with the ESMFold model running as an API
        in some remote server. 

        Parameters
        ----------
        url : str, optional
            The URL of the remote API to be called. The default is the URL
            for Meta's ESM Metagenomic Atlas (see [here](https://esmatlas.com/about#api)
            for more details).
        """
        self._url = url
        self._verify = self._url.find('esmatlas') == -1
        self._request_is_json = request_is_json
        self._response_is_json = response_is_json



    def predict_raw_pdb(self, sequence: str) -> str:
        """
        Calls the remote API at the given URL to predict the 3D structure
        of a given amino acid sequence.

        Parameters
        ----------
        sequence : str
            The amino acid sequence which structure will be predicted. Each residue
            must be represented with a single letter corresponding to one of the
            20 essential amino acids.
        
        Returns
        -------
        str
            The predicted 3D structure codified in a string as a PDB file.

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
        if self._url.find('esmatlas') != -1: time.sleep(1.5)
        if self._request_is_json:
            response = requests.post(self._url, 
                                    json={ 'input': sequence },
                                    timeout=30, 
                                    verify=self._verify)
        else:
            respose = requests.post(self._url, 
                                    data=sequence, 
                                    timeout=30, 
                                    verify=self._verify)
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
        if self._response_is_json:
            return response.json()['output']
        return response.content.decode()
