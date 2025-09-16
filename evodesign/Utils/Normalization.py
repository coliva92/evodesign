import numpy as np


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
