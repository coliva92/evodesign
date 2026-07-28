from ...System.APIRequester import APIRequester
from ...System.Exceptions import *
from .ESMFoldInterface import ESMFoldInterface




class ESMFoldRemoteAPI(ESMFoldInterface):

    def __init__(
        self,
        requester: APIRequester = APIRequester(
            url="https://api.esmatlas.com/foldSequence/v1/pdb/",
            # json_request_key="sequence",
            # json_response_key="pdb",
        ),
    ) -> None:
        super().__init__()
        self.requester = requester
        return



    def predict_single_pdb_str(
        self,
        sequence: str,
    ) -> str:
        return self.requester.post(sequence)
