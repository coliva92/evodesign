import numpy as np
import numpy.typing as npt


def reciprocal(
    x: float,
    scaling_factor: float = 1.0,
) -> float:
    return 1.0 / (1.0 + scaling_factor * x)


def sigmoid(
    x: float,
    scaling_factor: float = 1.0,
    offset: float = 0.0,
) -> float:
    return 1.0 / (1.0 + np.exp(scaling_factor * x + offset))


def z_score(
    x: float,
    mean_value: float,
    stdev: float,
) -> float:
    return (x - mean_value) / stdev


def cos_normalization(x: float, offset: float, power: float) -> float:
    return ((np.cos(x-offset)+1)/2)**power


def cos_similarity(u: npt.NDArray[np.float64], v: npt.NDArray[np.float64]) -> float:
    n1 = np.linalg.norm(u)
    n2 = np.linalg.norm(v)
    w = np.dot(u, v) / (n1*n2)
    return (1 + w) / 2
