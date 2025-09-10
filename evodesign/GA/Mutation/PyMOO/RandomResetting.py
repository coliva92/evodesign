from pymoo.core.mutation import Mutation
from pymoo.core.variable import Real
from pymoo.core.problem import Problem
from copy import deepcopy
import numpy as np
import numpy.typing as npt


class RandomResetting(Mutation):

    AMINO_ACIDS_ALPHABET = np.array(range(20), dtype=np.int64)

    def __init__(
        self,
        prob: float,
        prob_var: float,
    ) -> None:
        super().__init__()
        self.prob = Real(prob, bounds=(0.0, 1.0), strict=(0.0, 1.0))
        self.prob_var = Real(prob_var, bounds=(0.0, 1.0), strict=(0.0, 1.0))

    def _do(
        self,
        problem: Problem,
        X: npt.NDArray[np.int64],
        **kwargs,
    ):
        mask = np.random.random(X.shape) < self.prob_var.value
        mutations = np.random.choice(self.AMINO_ACIDS_ALPHABET[1:], size=X.shape)
        Xp = deepcopy(X)
        Xp[mask] = (X[mask] + mutations[mask]) % len(self.AMINO_ACIDS_ALPHABET)
        return Xp
