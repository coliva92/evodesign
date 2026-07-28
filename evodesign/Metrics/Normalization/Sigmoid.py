from .Formulas import sigmoid
from .Normalization import Normalization




class Sigmoid(Normalization):

    def __init__(
        self,
        scaling_factor: float = 1,
        offset: float = 0,
    ) -> None:
        super().__init__()
        self.scaling_factor = scaling_factor
        self.offset = offset
        return



    def do(self, x: float) -> float:
        return sigmoid(x, self.scaling_factor, self.offset)
