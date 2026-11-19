import numpy as np
import numpy.typing as npt
from scipy.special import rel_entr


def rmse(
    x: npt.NDArray[np.float64],
    y: npt.NDArray[np.float64],
) -> float:
    return np.sqrt(np.mean((x - y) ** 2))


def euclidean_distance(
    u: npt.NDArray[np.float64],
    v: npt.NDArray[np.float64],
) -> float:
    return np.linalg.norm(u - v)


def cosine_similarity(
    u: npt.NDArray[np.float64],
    v: npt.NDArray[np.float64],
) -> float:
    n1 = np.linalg.norm(u)
    n2 = np.linalg.norm(v)
    return np.dot(u, v) / (n1 * n2)


def norm_cosine_similarity(
    u: npt.NDArray[np.float64],
    v: npt.NDArray[np.float64],
) -> float:
    return (1 + cosine_similarity(u, v)) / 2.0


def jensen_shannon_distance(
    p: npt.NDArray[np.float64],
    q: npt.NDArray[np.float64],
) -> float:
    M = (p + q) / 2.0
    left = rel_entr(p, M)
    right = rel_entr(q, M)
    left_sum = np.max([0, np.sum(left)])
    right_sum = np.max([0, np.sum(right)])
    return np.sqrt((left_sum + right_sum) / 2.0)
