import numpy as np  # noqa: N999

from .Normalization import Normalization


class Sigmoid(Normalization):
    def __init__(
        self,
        scaling_factor: float = 1.0,
        offset: float = 0.0,
    ) -> None:
        super().__init__()
        self.scaling_factor = scaling_factor
        self.offset = offset

    def do(self, x: float) -> float:
        return 1.0 / (1.0 + np.exp(self.scaling_factor * x + self.offset))
