from pymoo.core.callback import Callback
from pymoo.core.algorithm import Algorithm as PyMOOAlgorithm
from ..Problems.MonoObjective.AlternatingMonoCPD import AlternatingMonoCPD
from typing import Optional


class MonoFitnessFnAlternator(Callback):

    def __init__(
        self,
        alt_every_nth_generation: int = 10,
        next_callback: Optional[Callback] = None,
    ):
        super().__init__()
        self.next_callback = next_callback
        self.alt_every_nth_generation = alt_every_nth_generation
        return

    def notify(self, algorithm: PyMOOAlgorithm) -> None:
        problem = algorithm.problem
        assert isinstance(problem, AlternatingMonoCPD)
        if (
            algorithm.n_gen % self.alt_every_nth_generation == 0
            or algorithm.n_gen % self.alt_every_nth_generation == 1
        ):
            problem.alternate_fitness_fn()
        if self.next_callback is not None:
            self.next_callback.notify(algorithm)
        return
