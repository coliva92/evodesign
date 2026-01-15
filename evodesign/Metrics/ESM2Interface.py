from ..RetrievableSettings import RetrievableSettings
from typing import Tuple
import numpy as np
import numpy.typing as npt


class ESM2Interface(RetrievableSettings):

    def query_model(
        self,
        sequence: str,
        sequence_name: str = "tmp_protein",
    ) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        raise NotImplementedError
