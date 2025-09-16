from .Normalization import Normalization
from ...Utils.Normalization import reciprocal


class Reciprocal(Normalization):

    def __init__(self, scaling_factor: float = 1):
        super().__init__()
        self.scaling_factor = scaling_factor

    def do(self, x: float) -> float:
        return reciprocal(x, self.scaling_factor)
