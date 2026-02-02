from pymoo.operators.sampling.rnd import IntegerRandomSampling
import numpy as np
import numpy.typing as npt


class ProfileIntegerSampling(IntegerRandomSampling):

    def _do(self, problem, n_samples, **kwargs):
        if problem.aa_profile is None:
            return super()._do(problem, n_samples)
        return self.generate_mutant_sequences(n_samples, problem.aa_profile)

    def generate_mutant_sequences(
        self,
        num_mutants: int,
        aa_profile: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.int64]:
        cumsum = aa_profile.cumsum(axis=1)
        # We reshape to (L, 1) to enable broadcasting against the (L, 20) cumsum
        rand_vals = np.random.rand(num_mutants, aa_profile.shape[0], 1)
        mutant_sequences = (cumsum <= rand_vals).sum(axis=-1)
        return mutant_sequences
