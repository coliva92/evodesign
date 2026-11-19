from .Normalization import Normalization  # noqa: N999


class ZScore(Normalization):
    def __init__(
        self,
        mean_value: float,
        stdev: float,
    ) -> None:
        super().__init__()
        self.mean_value = mean_value
        self.stdev = stdev

    def do(self, x: float) -> float:
        return (x - self.mean_value) / self.stdev
