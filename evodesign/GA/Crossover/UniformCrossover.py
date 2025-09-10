from .Crossover import Crossover
from pymoo.operators.crossover.ux import UX


class UniformCrossover(Crossover):

    def __init__(
        self,
        prob: float = 1.0,
    ):
        super().__init__(prob, UX(prob=prob))
