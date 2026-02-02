from ....Algorithms.ProfileIntegerSampling import ProfileIntegerSampling
from ....Utils.AminoAcids import AMINO_ACIDS_INT_ALPHABET
from pymoo.core.mutation import Mutation
from pymoo.core.variable import Real
from pymoo.core.problem import Problem
from copy import deepcopy
import numpy as np
import numpy.typing as npt


class RandomResetting(Mutation):

    def __init__(
        self,
        prob: float,
        prob_var: float,
    ) -> None:
        super().__init__()
        self.prob = Real(prob, bounds=(0.0, 1.0), strict=(0.0, 1.0))
        self.prob_var = Real(prob_var, bounds=(0.0, 1.0), strict=(0.0, 1.0))
        self.sampler = ProfileIntegerSampling()

    def _do(
        self,
        problem: Problem,
        X: npt.NDArray[np.int64],
        **kwargs,
    ):
        mask = np.random.random(X.shape) < self.prob_var.value
        Xp = deepcopy(X)
        if problem.aa_profile is None:
            mutations = np.random.choice(AMINO_ACIDS_INT_ALPHABET[1:], size=X.shape)
            Xp[mask] = (X[mask] + mutations[mask]) % len(AMINO_ACIDS_INT_ALPHABET)
        else:
            mutations = self.sampler.generate_mutant_sequences(
                X.shape[0], problem.aa_profile
            )
            Xp[mask] = mutations[mask]
        return Xp
