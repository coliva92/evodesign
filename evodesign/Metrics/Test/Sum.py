import numpy as np  # noqa: N999
import numpy.typing as npt

from ...Chemistry.Sequences import NUM_AMINO_ACIDS
from ..ContextInterface import ContextInterface
from ..SequenceBasedEvaluator import SequenceBasedEvaluator


class Sum(SequenceBasedEvaluator):
    """
    A dummy fitness function term that converts individual amino acid letters in a
    sequence into unique integer values and adds them together.
    `Sum` is meant to be used for making quick tests and check the code's integrity
    without needing to initialize costly deep learning models. As such, it does not
    represent true protein design.
    """

    _MAX_VALUE_PER_RESIDUE = NUM_AMINO_ACIDS - 1

    def do(
        self,
        model_sequence: npt.NDArray[np.int64],
        **kwargs,
    ) -> tuple[float, ...]:
        total_sum = model_sequence.sum(axis=0)
        avg = total_sum / (Sum._MAX_VALUE_PER_RESIDUE * model_sequence.shape[0])
        return total_sum, avg

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> dict[str, float]:
        model_sequence = context.get_model_chain().sequence_numpy
        total_sum, avg = self.do(model_sequence)
        return {"sum": total_sum, "average": avg}
