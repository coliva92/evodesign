from typing import List, Optional, Tuple

import numpy as np
import numpy.typing as npt

from ..RetrievableSettings import RetrievableSettings


class ESM2Interface(RetrievableSettings):

    def query_model(
        self,
        sequence: str,
        sequence_name: str,
        submap_indices: Optional[List[int]],
    ) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        raise NotImplementedError
