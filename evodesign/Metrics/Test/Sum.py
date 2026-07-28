from ..NonStructuralMetric import NonStructuralMetric
from ..ContextInterface import ContextInterface
from ...Chemistry.Sequences import NUM_AMINO_ACIDS
from typing import Dict, Tuple
import numpy as np
import numpy.typing as npt


class Sum(NonStructuralMetric):

    _MAX_VALUE_PER_RESIDUE = NUM_AMINO_ACIDS - 1

    def do(self, model_sequence: npt.NDArray[np.int64], **kwargs) -> Tuple[float, float]:
        total_sum = model_sequence.sum(axis=0)
        avg = total_sum / (Sum._MAX_VALUE_PER_RESIDUE * model_sequence.shape[0])
        return total_sum, avg

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> Dict[str, float]:
        model_sequence = context.get_model_chain().sequence_numpy
        total_sum, avg = self.do(model_sequence)
        return {
            "sum": total_sum,
            "average": avg
        }
