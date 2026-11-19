import numpy as np  # noqa: N999
import numpy.typing as npt

from ...System.APIRequester import APIRequester
from .ModelInterface import ModelInterface


class ESM2RemoteAPI(ModelInterface):
    def __init__(
        self,
        requester: APIRequester | None = None,
        **kwargs,
    ) -> None:
        super().__init__(**kwargs)
        self.requester = requester

    def query_model(
        self,
        sequence: str,
        sequence_name: str = "tmp_protein",
        submap_indices: list[int] | None = None,
    ) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:

        if self.requester is None:
            self.requester = APIRequester(
                url="http://127.0.0.1:8000/esm",
                json_request_key="sequence",
                json_response_key="results",
            )

        data = self.requester.post(sequence)
        molecular_descriptors = np.array(data["desc_matrix"])
        predicted_contact_map = np.array(data["predicted_contacts"])
        if submap_indices is not None:
            molecular_descriptors = molecular_descriptors[submap_indices]
            # convert the 1d array into a 2d array to simplify logic
            L = len(predicted_contact_map)
            # remember that L = N(N - 1) / 2, then apply the quadratic formula
            N = int((1 + np.sqrt(1 + 8 * L)) / 2)
            tmp_matrix = np.zeros((N, N))
            rows, cols = np.triu_indices(N, k=1)
            tmp_matrix[rows, cols] = predicted_contact_map
            tmp_matrix = tmp_matrix + tmp_matrix.T
            predicted_contact_map = tmp_matrix[np.ix_(submap_indices, submap_indices)]
        return molecular_descriptors, predicted_contact_map
