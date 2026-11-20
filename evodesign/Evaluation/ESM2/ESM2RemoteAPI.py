from __future__ import annotations  # noqa: N999

import numpy as np
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
        return molecular_descriptors, predicted_contact_map
