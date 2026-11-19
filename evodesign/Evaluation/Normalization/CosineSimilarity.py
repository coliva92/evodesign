from .Normalization import Normalization  # noqa: N999


class CosineSimilarity(Normalization):
    def do(self, x: float) -> float:
        return (1 + x) / 2.0
