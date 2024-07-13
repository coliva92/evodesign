from .Metric import Metric
from typing import Optional, List
from Bio.PDB.Atom import Atom
import pandas as pd
from ..Context import Context
import requests
from ..Exceptions import *
from .ESM2Descriptors import ESM2Descriptors
import time





class ESM2DescriptorsRemoteApi(Metric):

    def _params(self) -> dict:
        params = super()._params()
        params['url'] = self._url
        return params



    def __init__(self, 
                 url: str = "http://127.0.0.1:5000/esm",
                 sleep_time: float = 0.5,
                 column: Optional[str] = None
                 ) -> None:
        super().__init__(column)
        self._url = url
        self._sleep_time = sleep_time
    


    def _compute_values(self, 
                        backbone: List[Atom],
                        data: pd.Series,
                        context: Context
                        ) -> pd.Series:
        sequence, sequence_id = data["sequence"], data["sequence_id"]
        if self._sleep_time > 0.0:
            time.sleep(self._sleep_time)
        response = requests.post(self._url, 
                                 json={ "sequence": sequence },
                                 timeout=30)
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
        vectors = response.json()["descriptors"]
        txt_path = f"{context.workspace.esm2_dir}/{sequence_id}.txt"
        ESM2Descriptors.save_descriptor_vectors_txt(vectors, txt_path)
        data[self.column_name()] = txt_path
        return data
