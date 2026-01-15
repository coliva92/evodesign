from .ESM2Interface import ESM2Interface
from ..Utils.APIRequester import APIRequester
import numpy as np
import numpy.typing as npt
from typing import Tuple


class ESM2RemoteAPI(ESM2Interface):

    def __init__(
        self,
        requester: APIRequester = APIRequester(
            url="http://127.0.0.1:8000/esm",
            json_request_key="sequence",
            json_response_key="results",
        ),
    ) -> None:
        super().__init__()
        self.requester = requester
        return

    def query_model(
        self,
        sequence: str,
        sequence_name: str = "tmp_protein",
    ) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        data = self.requester.post(sequence)
        desc_matrix = np.array(data["desc_matrix"])
        predicted_contacts = np.array(data["predicted_contacts"])
        return desc_matrix, predicted_contacts
