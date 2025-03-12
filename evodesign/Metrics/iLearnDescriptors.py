from .Metric import Metric
import numpy.typing as npt
import numpy as np
from descproteins import *
from pubscripts import *
import numpy as np
import numpy.typing as npt


class iLearnDescriptors(Metric):

    def __init__(self, method: str) -> None:
        super().__init__()
        self.method = method

    def compute_descriptors_vector(
        self, fasta_path: str, sequence: str, sequence_name: str = "temp_protein"
    ) -> npt.NDArray[np.float64]:
        with open(fasta_path, "wt", encoding="utf-8") as fasta_file:
            fasta_file.write(f">1|{sequence_name}|{sequence_name}|\n{sequence}\n")
        sequences = read_fasta_sequences.read_fasta_sequences(fasta_path)
        encodings = eval(f"{self.method}.{self.method}(sequences)")
        return np.array(encodings[0])

    def do(
        self,
        model_fasta_path: str,
        model_sequence: str,
        ref_descriptors: npt.NDArray[np.float64],
        **kwargs,
    ) -> float:
        descriptors = self.compute_descriptors_vector(model_fasta_path, model_sequence)
        return np.sqrt(np.mean(descriptors - ref_descriptors) ** 2)
