from .Evaluator import Evaluator  # noqa: N999


class StructureBasedEvaluator(Evaluator):
    def requires_structure_predictor(self) -> bool:
        return True
