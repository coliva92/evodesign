from .Normalization import Normalization
import numpy as np


class Power(Normalization):

    def __init__(self, n_power: float):
        super().__init__()
        self.n_power = n_power
        return

    def do(self, x: float) -> float:
        return np.power(x, self.n_power)
