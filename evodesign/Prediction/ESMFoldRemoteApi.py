from .Predictor import Predictor
from typing import Optional
from ..Exceptions import *
import requests
import time


class ESMFoldRemoteApi(Predictor):

    def __init__(
        self,
        url: str = "https://api.esmatlas.com/foldSequence/v1/pdb/",
        json_request_key: Optional[str] = "sequence",
        json_response_key: Optional[str] = "pdb",
        sleep_time: float = 1.5,
        connection_timeout: int = 30,
        verify: bool = False,
    ) -> None:
        super().__init__()
        self.url = url
        self.json_request_key = json_request_key
        self.json_response_key = json_response_key
        self.sleep_time = sleep_time
        self.connection_timeout = connection_timeout
        self.verify = verify

    def predict_pdb_str(self, sequence: str) -> str:
        if self.sleep_time > 0.0:
            time.sleep(self.sleep_time)
        if self.json_request_key is not None:
            response = requests.post(
                self.url,
                json={self.json_request_key: sequence},
                timeout=self.connection_timeout,
                verify=self.verify,
            )
        else:
            response = requests.post(
                self.url,
                data=sequence,
                timeout=self.connection_timeout,
                verify=self.verify,
            )
        if response.status_code == 400:
            raise HttpBadRequest
        elif response.status_code == 403:
            raise HttpForbidden
        elif response.status_code == 504:
            raise HttpGatewayTimeout
        elif response.status_code == 500:
            raise HttpInternalServerError
        elif response.status_code != 200:
            print(response.status_code)
            print(response.content.decode())
            raise HttpUnknownError
        if self.json_response_key is not None:
            return response.json()[self.json_response_key]
        return response.content.decode()
