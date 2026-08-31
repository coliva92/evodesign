from .Normalization import Normalization
import numpy as np




class InvertedMinMax(Normalization):

    def __init__(self, min_val: float, max_val: float) -> None:
        super().__init__()
        self.min_val = min_val
        self.max_val = max_val
        self._denom = max_val - min_val
        return



    def do(self, x: float) -> float:
        norm = (self.max_val - x) / self._denom # move away from max_val
        result = max(0, min(1, norm)) # clamping
        return result
