import numpy as np
import numpy.typing as npt

from ..System.APIRequester import APIRequester
from .ESM2ModelInterface import ESM2ModelInterface


class ESM2RemoteAPI(ESM2ModelInterface):
    def __init__(
        self,
        requester: APIRequester | None = None,
    ) -> None:
        super().__init__()
        self.requester = requester
        if self.requester is None:
            self.requester = APIRequester(
                url="http://127.0.0.1:8000/esm",
                json_request_key="sequence",
                json_response_key="results",
            )

    def query_model(
        self,
        sequence: str,
        sequence_name: str = "tmp_protein",
        submap_indices: list[int] | None = None,
    ) -> tuple[npt.NDArray[np.float64], ...]:
        data = self.requester.post(sequence)
        desc_matrix = np.array(data["desc_matrix"])
        predicted_contacts = np.array(data["predicted_contacts"])
        if submap_indices is not None:
            desc_matrix = desc_matrix[submap_indices]
            # convert the 1d array into a 2d array to simplify logic
            L = len(predicted_contacts)
            N = int((1 + np.sqrt(1 + 8 * L)) / 2)
            tmp_matrix = np.zeros((N, N))
            rows, cols = np.triu_indices(N, k=1)
            tmp_matrix[rows, cols] = predicted_contacts
            tmp_matrix = tmp_matrix + tmp_matrix.T
            predicted_contacts = tmp_matrix[np.ix_(submap_indices, submap_indices)]
        return desc_matrix, predicted_contacts
