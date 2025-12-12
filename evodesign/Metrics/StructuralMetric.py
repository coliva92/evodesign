from .Metric import Metric


class StructuralMetric(Metric):

    def requires_structure_predictor(self) -> bool:
        return True
