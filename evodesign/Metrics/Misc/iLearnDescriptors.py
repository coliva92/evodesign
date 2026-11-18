import numpy as np  # noqa: N999
import numpy.typing as npt
from descproteins import *  # assuming iLearn source code is globally accessible
from pubscripts import *

from ..ContextInterface import ContextInterface
from ..metrics import rmse
from ..NonStructuralMetric import NonStructuralMetric


class iLearnDescriptors(NonStructuralMetric):
    def __init__(
        self,
        method: str,
    ) -> None:
        super().__init__()
        self.method = method

    def compute_descriptors_vector(
        self,
        fasta_filename: str,
        sequence: str,
        sequence_name: str = "tmp_protein",
    ) -> npt.NDArray[np.float64]:
        with open(fasta_filename, "wt", encoding="utf-8") as fasta_file:
            fasta_file.write(f">1|{sequence_name}|{sequence_name}|\n{sequence}\n")
        _ = read_fasta_sequences.read_fasta_sequences(fasta_filename)
        encodings = eval(f"{self.method}.{self.method}(sequences)")
        return np.array(encodings[0])

    def do(
        self,
        model_descriptors: npt.NDArray[np.float64],
        ref_descriptors: npt.NDArray[np.float64],
        **kwargs,
    ) -> float:
        descriptors = self.compute_descriptors_vector(model_fasta_path, model_sequence)
        return rmse(descriptors - ref_descriptors)

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> dict[str, float]:
        fasta_filename = "tmp_sequence.fasta"
        ref_descriptors = context.get_extra_param_value(
            f"reference_ilearn_desc_{self.method}"
        )
        if ref_descriptors is None:
            ref_sequence = context.get_reference_chain().sequence
            ref_descriptors = self.compute_descriptors_vector(
                fasta_filename,
                ref_sequence,
            )
            context.set_extra_param_value(
                ref_descriptors,
                f"reference_ilearn_desc_{self.method}",
            )
        model_sequence = context.get_model_chain().sequence
        model_descriptors = self.compute_descriptors_vector(
            fasta_filename,
            model_sequence,
        )
        rmse = self.do(model_descriptors, ref_descriptors)
        return {"rmse": rmse}
