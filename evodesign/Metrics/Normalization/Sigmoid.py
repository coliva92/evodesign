from .Normalization import Normalization
from ...Utils.Normalization import sigmoid


class Sigmoid(Normalization):

    def __init__(
        self,
        scaling_factor: float = 1,
        offset: float = 0,
    ):
        super().__init__()
        self.scaling_factor = scaling_factor
        self.offset = offset

    def do(self, x: float) -> float:
        return sigmoid(x, self.scaling_factor, self.offset)
