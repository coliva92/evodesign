from .Metric import Metric
from typing import Optional, List
from Bio.PDB.Atom import Atom
import pandas as pd
from ..Context import Context
import requests
from ..Exceptions import *
from .ESM2Descriptors import ESM2Descriptors
import time
import numpy as np
import numpy.typing as npt





class ESM2DescriptorsRemoteApi(ESM2Descriptors):

    def _params(self) -> dict:
        return {
            "column": self._column_name,
            "url": self._url,
            "sleep_time": self._sleep_time
        }



    def __init__(self, 
                 url: str = "http://127.0.0.1:5000/esm",
                 sleep_time: float = 0.5,
                 column: Optional[str] = None
                 ) -> None:
        super().__init__(column)
        self._url = url
        self._sleep_time = sleep_time
    


    def compute_descriptor_vectors(self, 
                                   sequence: str,
                                   sequence_id: str = 'prot_0000_0000'
                                   ) -> npt.NDArray[np.float64]:
        if self._sleep_time > 0.0:
            time.sleep(self._sleep_time)
        response = requests.post(self._url, 
                                 json={ "sequence": sequence },
                                 timeout=120,
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
        return response.json()["descriptors"]
