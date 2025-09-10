from .Metric import Metric
from .ContextInterface import ContextInterface
from typing import Dict
import numpy.typing as npt
import numpy as np
from descproteins import *  # assuming iLearn source code is globally accessible
from pubscripts import *
import numpy as np
import numpy.typing as npt
import os


class iLearnDescriptorsRMSE(Metric):

    def __init__(self, method: str) -> None:
        super().__init__()
        self.method = method

    def do(
        self,
        model_fasta_path: str,
        model_sequence: str,
        ref_descriptors: npt.NDArray[np.float64],
        **kwargs,
    ) -> float:
        descriptors = self.compute_descriptors_vector(model_fasta_path, model_sequence)
        return np.sqrt(np.mean(descriptors - ref_descriptors) ** 2)

    def do_for_fitness_fn(
        self,
        context: ContextInterface,
    ) -> Dict[str, float]:
        fasta_path = context.get_extra_param_value("fasta_path")
        if fasta_path is None:
            workspace_dir = context.get_extra_param_value("workspace_dir")
            if workspace_dir is None:
                raise KeyError
            fasta_path = os.path.join(workspace_dir, "tmp_sequence.fasta")
            context.set_extra_param_value("fasta_path", fasta_path)
        ref_descriptors = context.get_extra_param_value(
            f"reference_ilearn_desc_{self.method}"
        )
        if ref_descriptors is None:
            ref_sequence = context.get_reference_chain().sequence
            ref_descriptors = self.compute_descriptors_vector(fasta_path, ref_sequence)
            context.set_extra_param_value(
                ref_descriptors, f"reference_ilearn_desc_{self.method}"
            )
        model_sequence = context.get_model_chain().sequence
        rmse = self.do(fasta_path, model_sequence, ref_descriptors)
        return {"rmse": rmse}

    def compute_descriptors_vector(
        self,
        fasta_path: str,
        sequence: str,
        sequence_name: str = "tmp_protein",
    ) -> npt.NDArray[np.float64]:
        with open(fasta_path, "wt", encoding="utf-8") as fasta_file:
            fasta_file.write(f">1|{sequence_name}|{sequence_name}|\n{sequence}\n")
        sequences = read_fasta_sequences.read_fasta_sequences(fasta_path)
        encodings = eval(f"{self.method}.{self.method}(sequences)")
        return np.array(encodings[0])
