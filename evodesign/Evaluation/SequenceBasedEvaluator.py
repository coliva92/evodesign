from .Evaluator import Evaluator  # noqa: N999


class SequenceBasedEvaluator(Evaluator):
    def requires_structure_prediction(self) -> bool:
        return False
