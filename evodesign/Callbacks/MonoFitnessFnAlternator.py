from pymoo.core.algorithm import Algorithm as PyMOOAlgorithm
from pymoo.core.callback import Callback

from ..Problems.MonoObjective.AlternatingMonoCPD import AlternatingMonoCPD




class MonoFitnessFnAlternator(Callback):

    def __init__(
        self,
        alt_every_nth_gen: int = 10,
    ):
        super().__init__()
        self.alt_every_nth_gen = alt_every_nth_gen
        return



    def notify(self, algorithm: PyMOOAlgorithm) -> None:
        problem = algorithm.problem
        assert isinstance(problem, AlternatingMonoCPD)
        if (
            algorithm.n_gen % self.alt_every_nth_gen == 0
            or algorithm.n_gen % self.alt_every_nth_gen == 1
        ):
            problem.alternate_fitness_fn()
        return
