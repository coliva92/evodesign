from .Metric import Metric


class NonStructuralMetric(Metric):

    def requires_structure_predictor(self) -> bool:
        return False
