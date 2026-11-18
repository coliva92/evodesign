import numpy as np  # noqa: N999

from .Normalization import Normalization


class Gaussian(Normalization):
    def __init__(
        self,
        mean_val: float = 1.0,
        stdev: float = 1.0,
        scaling_factor: float = 1.0,
    ) -> None:
        super().__init__()
        self.mean_val = mean_val
        self.stdev = stdev

    def do(self, x: float) -> float:
        return np.exp(
            -((x - self.mean_val) ** 2) / (2.0 * self.scaling_factor * self.stdev**2)
        )
