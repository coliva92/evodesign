from .Formulas import reciprocal
from .Normalization import Normalization




class Reciprocal(Normalization):

    def __init__(self, scaling_factor: float = 1) -> None:
        super().__init__()
        self.scaling_factor = scaling_factor
        return



    def do(self, x: float) -> float:
        return reciprocal(x, self.scaling_factor)
