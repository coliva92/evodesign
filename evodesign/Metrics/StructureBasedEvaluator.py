from .Metric import Metric  # noqa: N999


class StructureBasedEvaluator(Metric):
    def requires_structure_predictor(self) -> bool:
        return True
