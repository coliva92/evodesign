from .Normalization import Normalization


class MinMax(Normalization):

    def __init__(self, min_val: float, max_val: float) -> None:
        super().__init__()
        self.min_val = min_val
        self.max_val = max_val
        self._denom = max_val - min_val
        return

    def do(self, x: float) -> float:
        return (x - self.min_val) / self._denom
