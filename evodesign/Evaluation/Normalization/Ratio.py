from .Normalization import Normalization  # noqa: N999


class Ratio(Normalization):
    def __init__(self, scaling_factor: float = 1.0):
        super().__init__()
        self.scaling_factor = scaling_factor

    def do(self, x: float) -> float:
        return 1.0 - (x / self.scaling_factor)
