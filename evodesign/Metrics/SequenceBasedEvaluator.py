from .Metric import Metric  # noqa: N999


class SequenceBasedEvaluator(Metric):
    def requires_structure_predictor(self) -> bool:
        return False
