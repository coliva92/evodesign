from .Normalization import Normalization  # noqa: N999


class Reciprocal(Normalization):
    def __init__(
        self,
        scaling_factor: float = 1.0,
    ) -> None:
        super().__init__()
        self.scaling_factor = scaling_factor

    def do(self, x: float) -> float:
        return 1.0 / (1.0 + self.scaling_factor * x)
