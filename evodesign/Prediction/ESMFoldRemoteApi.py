from .Predictor import Predictor
from ..Utils.APIRequester import APIRequester
from ..Exceptions import *


class ESMFoldRemoteApi(Predictor):

    def __init__(
        self,
        requester: APIRequester = APIRequester(
            url="https://api.esmatlas.com/foldSequence/v1/pdb/",
            json_request_key="sequence",
            json_response_key="pdb",
        ),
    ) -> None:
        super().__init__()
        self.requester = requester

    def predict_pdb_str(self, sequence: str) -> str:
        return self.requester.post(sequence)
