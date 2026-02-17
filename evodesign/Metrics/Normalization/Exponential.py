from .Normalization import Normalization
import numpy as np


class Exponential(Normalization):

    def __init__(self, scaling_factor: float = 1.0) -> None:
        super().__init__()
        self.scaling_factor = scaling_factor
        return

    def do(self, x: float) -> float:
        return (np.exp(self.scaling_factor * x) - 1) / (np.exp(self.scaling_factor) - 1)
