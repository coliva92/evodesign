from .Evaluator import Evaluator  # noqa: N999


class StructureBasedEvaluator(Evaluator):
    def requires_structure_prediction(self) -> bool:
        return True
