from .Predictor import Predictor
from ..Exceptions import *
import requests
import time





class ESMFoldRemoteApi(Predictor):

    def _params(self) -> dict:
        parent_params = super()._params()
        params = { 
            'url': self._url,
            'request_is_json': self._request_is_json,
            'response_is_json': self._response_is_json,
            'sequence_json_key': self._sequence_json_key,
            'prediction_json_key': self._prediction_json_key,
            'sleep_time': self._sleep_time
        }
        return { **parent_params, **params }
    


    def __init__(self, 
                 url: str = 'https://api.esmatlas.com/foldSequence/v1/pdb/',
                 request_is_json: bool = False,
                 response_is_json: bool = False,
                 sequence_json_key: str = 'sequence',
                 prediction_json_key: str = 'pdb',
                 sleep_time: float = 1.5
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
        request_is_json : bool, optional
            Indicates if the request data is sent to the remote API in JSON 
            format or as a raw string in the request body. Default is `False`.
        response_is_json : bool, optional
            Indicates if the response data is recieved from the remote API in 
            JSON format or as a raw string in the request body. Default is 
            `False`.
        sequence_json_key : str, optional
            If `request_is_json` is `True`, then `sequence_json_key` describes 
            the name of the field in the request JSON that will contain the 
            amino acid sequence. Default is "input".
        prediction_json_key : str, optional
            If `response_is_json` is `True`, then `prediction_json_key` 
            describes the name of the field in the response JSON that should
            contain the PDB string of the predicted structure. Default is 
            "output".
        sleep_time : float, optional
            A waiting time (in seconds) set in between requests to the remote 
            API. Increase this number if the remote API limites the rate at 
            which it admits requests. Set to zero to send requests with no wait
            time between them. Default is 1.5 seconds.
        """
        super().__init__()
        self._url = url
        self._request_is_json = request_is_json
        self._response_is_json = response_is_json
        self._sequence_json_key = sequence_json_key
        self._prediction_json_key = prediction_json_key
        self._sleep_time = sleep_time
        self._verify = self._url.find('esmatlas') == -1



    def predict_pdb_str(self, sequence: str) -> str:
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
        if self._sleep_time > 0: 
            time.sleep(self._sleep_time)
        if self._request_is_json:
            response = requests.post(self._url, 
                                    json={ self._sequence_json_key: sequence },
                                    timeout=30, 
                                    verify=self._verify)
        else:
            response = requests.post(self._url, 
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
            return response.json()[self._prediction_json_key]
        return response.content.decode()
