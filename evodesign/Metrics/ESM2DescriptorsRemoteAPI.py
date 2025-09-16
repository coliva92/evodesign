from .ESM2Descriptors import ESM2Descriptors
from ..Utils.APIRequester import APIRequester
import numpy as np
import numpy.typing as npt
from .Normalization.Normalization import Normalization
from .Normalization.Reciprocal import Reciprocal
from typing import Optional


class ESM2DescriptorsRemoteAPI(ESM2Descriptors):

    def __init__(
        self,
        requester: APIRequester = APIRequester(
            url="http://127.0.0.1:8000/esm",
            json_request_key="sequence",
            json_response_key="descriptors",
        ),
        normalization: Optional[Normalization] = Reciprocal(),
    ) -> None:
        super().__init__(
            normalization=normalization,
            gpu_device=None,
        )
        self.requester = requester

    def compute_descriptors_matrix(
        self,
        sequence: str,
        sequence_name: str = "tmp_protein",
    ) -> npt.NDArray[np.float64]:
        return np.array(self.requester.post(sequence))
