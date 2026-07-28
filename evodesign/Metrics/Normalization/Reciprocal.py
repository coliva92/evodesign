from .Formulas import reciprocal
from .Normalization import Normalization




class Reciprocal(Normalization):

    def __init__(self, scaling_factor: float = 1):
        super().__init__()
        self.scaling_factor = scaling_factor



    def do(self, x: float) -> float:
        return reciprocal(x, self.scaling_factor)
