import numpy as np
import numpy.typing as npt


AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
AMINO_ACIDS_INT_ALPHABET = np.array(range(20), dtype=np.int64)
MAP_AMINO_ACID_TO_INT = {aa: i for i, aa in enumerate(AMINO_ACIDS)}


def to_numpy(sequence: str) -> npt.NDArray[np.int64]:
    return np.array(
        [MAP_AMINO_ACID_TO_INT[aa] for aa in sequence], dtype=np.int64
    )


def to_str(sequence_numpy: npt.NDArray[np.int64]) -> str:
    return "".join(AMINO_ACIDS[i] for i in sequence_numpy)
