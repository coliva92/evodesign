from .Evaluator import Evaluator  # noqa: N999


class SequenceBasedEvaluator(Evaluator):
    def requires_structure_predictor(self) -> bool:
        return False
