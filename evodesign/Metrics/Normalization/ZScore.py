from .Formulas import z_score
from .Normalization import Normalization


class ZScore(Normalization):

    def __init__(
        self,
        mean_value: float,
        stdev: float,
    ):
        super().__init__()
        self.mean_value = mean_value
        self.stdev = stdev

    def do(self, x: float) -> float:
        return z_score(x, self.mean_value, self.stdev)
