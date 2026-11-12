import numpy as np
import numpy.typing as npt


class ESM2Interface:
    def query_model(
        self,
        sequence: str,
        sequence_name: str,
        submap_indices: list[int] | None,
    ) -> tuple[npt.NDArray[np.float64], ...]:
        raise NotImplementedError
